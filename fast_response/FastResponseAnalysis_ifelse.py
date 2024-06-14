r''' General Fast Response Analysis Class.

    Author: Alex Pizzuto
    Date: 2021

    One of two variants to include mutiple datasets
    '''

from abc import abstractmethod
import os, sys, time, subprocess
import pickle, dateutil.parser, logging, warnings
from argparse import Namespace

import h5py
import healpy                 as hp
import numpy                  as np
import seaborn                as sns
import matplotlib             as mpl
import matplotlib.pyplot      as plt
import numpy.lib.recfunctions as rf
from astropy.time             import Time
from scipy.special            import erfinv
from matplotlib.lines         import Line2D

from skylab.datasets          import Datasets
from skylab.llh_models        import EnergyLLH
from skylab.priors            import SpatialPrior
from skylab.ps_injector       import PointSourceInjector
from skylab.ps_llh            import PointSourceLLH, MultiPointSourceLLH
from skylab.ps_injector       import PriorInjector
from skylab.spectral_models   import PowerLaw 
from skylab.temporal_models   import BoxProfile, TemporalModel
# import meander

from . import web_utils
from . import sensitivity_utils
from . import plotting_utils
from .reports import FastResponseReport

mpl.use('agg')
current_palette = sns.color_palette('colorblind', 10)
logging.getLogger().setLevel(logging.ERROR)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

class FastResponseAnalysis(object):
    """ 
    Object to do realtime followup analyses of 
    astrophysical transients with arbitrary event
    localization
    """
    _dataset = None
    _fix_index = True
    _float_index = not _fix_index
    _index_range = [1., 4.]
    _index = 2.0
    _floor = np.radians(0.2)
    _verbose = True
    _angScale = 2.145966
    _llh_seed = 1
    _season_names = [f"IC86, 20{y:02d}" for y in range(11, 22+1)]
    _nb_days = 10
    _ncpu = 5

    def __init__(self, name, tstart, tstop, skipped=None, seed=None,
                 outdir=None, save=True, extension=None):
        self.name = name
        self.multi = type(self.dataset) is list
        
        if seed is not None:
            self.llh_seed(seed)
        if outdir is None:
            outdir = os.environ.get('FAST_RESPONSE_OUTPUT')
            if outdir is None:
                outdir = os.getcwd()
        self.outdir = outdir

        start = Time(dateutil.parser.parse(tstart)).mjd
        stop = Time(dateutil.parser.parse(tstop)).mjd

        start_date = Time(dateutil.parser.parse(tstart)).datetime
        start_str = f'{start_date.year:02d}_' \
            + f'{start_date.month:02d}_{start_date.day:02d}'
        
        self.start = start
        self.stop = stop
        self.duration = stop - start
        self.centertime = (start + stop) / 2.
        
        self.analysisid = start_str + '_' + self.name.replace(' ', '_') 
        self.analysispath = os.path.join(self.outdir, self.analysisid)
        self.save_output = save

        if os.path.isdir(self.analysispath) and self.save_output:
            print("Directory {} already exists. Deleting it ...".format(
                self.analysispath))
            subprocess.call(['rm', '-r', self.analysispath])
            subprocess.call(['mkdir', self.analysispath])
            # sys.exit()
        elif self.save_output:
            subprocess.call(['mkdir', self.analysispath])

        if 'test' in self.name.lower():
            self.scramble = True
            if self._verbose:
                print('Working on scrambled data')
        else:
            self.scramble = False
            if self._verbose:
                print('Working on unscrambled (UNBLINDED) data')

        self.save_items = {
            'start':self.start, 'stop': self.stop,
            'name': self.name, 'analysisid': self.analysisid,
            'analysispath': self.analysispath
        }
        self.extension = extension
        if self.extension is not None:
            self.extension = np.deg2rad(self.extension)
        self.llh_seed = seed if seed is not None else 1
        self.skipped = skipped
        self.skipped_event = None
        self.exp = None
        self.spectrum = None
        if self._fix_index:
            self.spectrum = PowerLaw(A=1, gamma=self.index, E0=1000.)
        self.time_profile = BoxProfile(self.start, self.stop)
        self.llh = self.initialize_llh(skipped=skipped, scramble=self.scramble)
        self.inj = None

    @property
    def dataset(self):
        """Returns the datasets used"""
        return self._dataset
    @dataset.setter
    def dataset(self, x):
        self._dataset = x

    @property
    def index(self):
        """Returns the spectral index"""
        return self._index
    @index.setter
    def index(self, x):
        self._index = x

    @property
    def llh_seed(self):
        """Returns the seed used in the LLH"""
        return self._llh_seed
    @llh_seed.setter
    def llh_seed(self, x):
        self._llh_seed = x

    @property
    def llh_exp(self):
        """Returns a flat array of experimental data loaded across the used dataset(s) loaded into the LLH,
        limited to fields used for plotting and amended with a field for the dataset."""
        # NOTE can not use self.dset_container as that is before self.remove_event(skipped=skipped)
        # TODO make this more efficient?
        if self.multi:
            exp = {}
            for enum, _llh in self.llh._samples.items():
                exp[enum] = _llh.exp
        else:
            exp = {-1: self.llh.exp}
            
        merged_dtype = np.dtype([('run', '<i4'), ('event', '<i4'), ('time', '<f4'), ('ra', '<f4'), ('dec', '<f4'), ('sigma', '<f4'), ('enum', '<i4')])
        for enum, _exp in exp.items():
            _exp = rf.drop_fields(_exp, [_field for _field in _exp.dtype.names if _field not in merged_dtype.names])
            _exp = rf.append_fields(_exp, 'enum', np.full( _exp.size, enum))
            _exp = _exp.astype(merged_dtype)
            exp[enum] = _exp
        return np.concatenate([exp[enum] for enum in exp])



    def get_data_single(self, dataset, livestream_start=None, livestream_stop=None):
        """
        Gets the skylab data and MC from querying the i3live livestream for a single dataset.

        Parameters
        -----------
        dataset: str
            Name of dataset in skylab.datasets.Datasets
        livestream_start: float
            (optional) start time for getting data (MJD) from i3Live. Default is start time of the analysis - 5 days
        livestream_stop: float
            (optional) stop time for getting data (MJD). Default is stop time of the analysis

        Returns:
        --------
        dset_container: argparse.Namespace
            Namespace containing
            exp, mc, grl = single arrays of experimental data, Monte Carlo, GRL
            livetime = sum of livetime
            dset = Skylab dataset object
            sinDec_bins, energy_bins = binning to be used for LLH
        """
        if self._verbose:
            print(f"Grabbing data for {dataset}")



        dset = Datasets[dataset]
        dset_container = Namespace()
        # not all datasets have seasons from 2011--2022
        # e.g. the default GFU is 2011-2019, Greco 2012--2022
        archival_seasons = sorted(list(set(dset.season_names()).intersection(self._season_names)))
        
        # TODO change this if the used GFU ever gets updated, or replace with and of GRL
        if self.stop < 59215:
            if self._verbose:
                print("Old times, just grabbing archival data")
            exps, grls = [], []
            for season in archival_seasons:
                exp, mc, livetime = dset.season(season, floor=self._floor)
                grl = dset.grl(season)
                exps.append(exp)
                grls.append(grl)
            # TODO this relies on the assumption the archival GFU is equal to the default
            if (self.stop > 58933.0) and dataset.startswith('GFUOnline'):
                if self._verbose:
                    print('adding 2020 GFU from /data/user/apizzuto')
                # Add local 2020 if need be
                # TODO: Need to figure out what to do for zenith_smoothed
                exp_new = np.load(
                    '/data/user/apizzuto/fast_response_skylab/fast-response/'
                    + 'fast_response/2020_data/IC86_2020_data.npy')
                grl = np.load(
                    '/data/user/apizzuto/fast_response_skylab/fast-response/' 
                    + 'fast_response/2020_data/GRL/IC86_2020_data.npy')
                exp_new.dtype.names = exp.dtype.names
                exp_new['sigma'] = np.where(
                    exp_new['sigma'] > self._floor, 
                    exp_new['sigma'], 
                    self._floor)
                exps.append(exp_new)
                grls.append(grl)
            exp = np.concatenate(exps)
            grl = np.concatenate(grls)
        # TODO add livestream season to Greco dataset - big task
        else:
            if self._verbose:
                print("Recent time: querying the i3live database")
            if livestream_start is None or livestream_stop is None:
                livestream_start = self.start - 6.
                livestream_stop = self.stop
            exp, mc, livetime, grl = dset.livestream(
                livestream_start, livestream_stop,
                append=archival_seasons, 
                floor=self._floor)
        exp.sort(order='time')
        grl.sort(order='run')
        livetime = grl['livetime'].sum()
        # workaround while not all datasets have livestream yet
        if 'livestream' in dset.season_names():
            reference_season = 'livestream'
        else:
            reference_season = archival_seasons[0]
        sinDec_bins = dset.sinDec_bins(reference_season)
        energy_bins = dset.energy_bins(reference_season)

        dset_container.exp = exp
        dset_container.mc = mc
        dset_container.grl = grl
        dset_container.livetime = livetime
        dset_container.dset = dset
        dset_container.sinDec_bins = sinDec_bins
        dset_container.energy_bins = energy_bins
        return dset_container

    def get_data(self, **kwargs):
        """
        Gets the skylab data and MC from querying the i3live livestream

        Parameters
        -----------
        livestream_start: float
            (optional) start time for getting data (MJD) from i3Live. Default is start time of the analysis - 5 days
        livestream_stop: float
            (optional) stop time for getting data (MJD). Default is stop time of the analysis
        """
        if self.multi:
            dset_multi = {enum:self.get_data_single(_dataset) for enum, _dataset in enumerate(self.dataset)}
            self.dset_container = dset_multi
            # pivot from dictionary of namespaces
            # to dictionaries in the class namespace
            for enum, dset in dset_multi.items():
                for attr in vars(dset):
                    if not hasattr(self, attr) or getattr(self, attr) is None:
                        setattr(self, attr, {})
                    getattr(self, attr)[enum] = getattr(dset, attr)
        else:
            dset = self.get_data_single(self.dataset, **kwargs)
            self.dset_container = {-1: dset}
            for attr in vars(dset):
                setattr(self, attr, getattr(dset, attr))

        

        
    def initialize_llh(self, skipped=None, scramble=False):
        """
        Grab data and format it all into a skylab llh object

        Parameters
        -----------
        skipped: array of tuples 
            (optional) event(s) to be removed in the analysis. 
            Format: [(run_id,event_id),...]
        scramble: bool
            (optional) run on scrambled data (default=False)

        Returns
        ----------
        llh: skylab llh object
            Likelihood for the analysis
        """
        if self.exp is None:
            self.get_data()

        if self._verbose:
            print("Initializing Point Source LLH in Skylab")

        assert self._fix_index != self._float_index,\
            'Must choose to either float or fix the index'

        llh_list = []
        # kwargs to be broadcast to the BaseLLH constructor
        # not specific to dataset
        base_kwargs = dict(
            nsource=1.,                    # seed for nsignal fit
            nsource_bounds=(0., 1e3),      # bounds on fitted ns
            ncpu=self._ncpu,               # use 10 CPUs when computing trials
            seed=self.llh_seed,
        )

        for enum in self.dset_container:
            dset = self.dset_container[enum]
            twodim_bins = [dset.energy_bins, dset.sinDec_bins]
            if self._fix_index:
                llh_model = EnergyLLH(
                    twodim_bins=twodim_bins,
                    allow_empty=True,
                    spectrum=self.spectrum,
                    ncpu=self._ncpu) 
            elif self._float_index:
                llh_model = EnergyLLH(
                    twodim_bins=twodim_bins,
                    allow_empty=True,
                    bounds=self._index_range,
                    seed = self.index,
                    ncpu=self._ncpu)
            
            box = TemporalModel(
                grl=dset.grl,
                poisson_llh=True,
                days=self._nb_days,
                signal=self.time_profile)
            
            if skipped is not None:
                dset.exp = self.remove_event(dset.exp, dset.dset, skipped)

            if self.extension is not None:
                src_extension = float(self.extension)
                self.save_items['extension'] = src_extension
            else:
                src_extension = None

            llh = PointSourceLLH(
                dset.exp,                      # array with data events
                dset.mc,                       # array with Monte Carlo events
                dset.livetime,                 # total livetime of the data events
                scramble=scramble,             # set to False for unblinding
                timescramble=True,             # not just RA scrambling
                llh_model=llh_model,           # likelihood model
                temporal_model=box,            # use box for temporal model
                src_extension = src_extension, # Model symmetric extended source
                **base_kwargs,
                )
            llh_list.append(llh)

        if self.multi:
            multi_llh = MultiPointSourceLLH(**base_kwargs)
            for name, llh in zip(self.dataset, llh_list):
                llh.do_trials_seed = 1
                multi_llh.add_sample(name, llh)
            return multi_llh
        
        elif len(llh_list)==0:
            return llh_list[0]
        
        else:
            raise AssertionError('have multiple instances of PointSourceLLH but FRA is not multi')
    
    def remove_event(self, exp, dset, skipped):
        """
        Remove a given event from the analysis, eg. for an alert event

        Parameters
        -----------
        exp: skylab data object
            Data events loaded from livestream
        dset: skylab Dataset object
            Dataset used in the analysis
        skipped: array of tuples 
            Event(s) to be attempted to be removed in the analysis. 
            Format: [(run_id,event_id),...]

        Returns
        ----------
        exp: skylab data object
            Data events, minus the removed event (if successful)
        """
        try:
            event = skipped[0]
            mjd_keys = exp['time'][(exp['run'] == int(event[0])) & (exp['event'] == int(event[1]))]
            self.skipped_event = exp[exp['time'] == mjd_keys[0]][0]
            exp = dset.remove_ev(exp, mjd_keys=mjd_keys[0])
            self.save_items['skipped_event'] = self.skipped_event
        except:
            print("Was not able to remove event {}".format(skipped))
        return exp

    @abstractmethod
    def initialize_injector(self):
        pass

    @abstractmethod
    def unblind_TS(self):
        pass

    @abstractmethod
    def run_background_trials(self, ntrials=1000):
        pass

    def calc_pvalue(self, ntrials=1000, run_anyway=False):
        r""" Given an unblinded TS value, calculate the p value

        Parameters
        -----------
        ntrials: int
            Number of trials to run (default 1000)
        run_anyway: bool
            Choose to override and run background trials, even if TS=0 (default False)

        Returns
        --------
        p: float
            P value for this analysis
        """ 
        if self.ts is None:
            if self._verbose:
                print("Need TS value to find p value")
            self.unblind_TS()
        if self.ts == 0.0 and not run_anyway:
            if self._verbose:
                print("TS=0, no need to run background trials")
            p = 1.0
            self.tsd = None
        else:
            self.run_background_trials(ntrials=ntrials)
            p = np.count_nonzero(self.tsd >= self.ts) / float(self.tsd.size)
            sigma = self.significance(p)
            if self._verbose:
                print(f"p: {p:.4f}")
                print("Significance: {:.3f}sigma\n\n".format(sigma))
        self.p = p
        self.save_items['p'] = p
        return p

    @abstractmethod
    def find_coincident_events(self):
        pass

    @abstractmethod
    def upper_limit(self):
        pass

    @abstractmethod
    def write_circular(self):
        pass

    @abstractmethod
    def make_dNdE(self):
        pass

    @abstractmethod
    def make_all_report_plots(self):
        pass

    @abstractmethod
    def make_all_report_tables(self):
        pass

    def plot_tsd(self, allow_neg=False):
        r"""
        Outputs a plot of the background TS distributions
        as well as the observed TS

        Parameters
        -----------
        allow_neg: bool
            Choose to allow negative TS values (all negative values are then TS=0). Default False
        """
        if self.tsd is not None:
            fig, ax = plt.subplots()
            if allow_neg:
                if str(np.min(self.tsd))== '-inf': 
                    lower=-500.
                else: 
                    lower = np.max([np.min(self.tsd), -500.])

                #lowest bin as underflow bin
                self.tsd[self.tsd < lower] = lower 
                self.tsd[np.isinf(self.tsd)]= lower

                pos_bins=np.linspace(0., 25., 30) #fine pos bins
                neg_bins=np.linspace(lower, 0., 30) #coarse neg bins

                bins=np.concatenate([ neg_bins[:-1], pos_bins ]) 
            else:
                if (np.min(self.tsd) < 0.0) or (str(np.min(self.tsd))== '-inf'):
                    self.tsd[self.tsd < 0.] = 0.
                    self.tsd[np.isinf(self.tsd)]= 0.
                bins = np.linspace(0., 25., 30)
            
            plt.hist(self.tsd, bins= bins, 
                    label="Background Scrambles", density=True)
            if self.ts >= -500.:
                plt.axvline(self.ts, color = 'k', label = "Observed TS")
            else: 
                plt.annotate("Unblinded TS: {:.0f}".format(self.ts), (lower+20, 1e-3))
            plt.legend(loc=1)

            if bins[0]<0.:
                plt.xlabel("TS (Note: lowest bin is an underflow bin)")
            else:
                plt.xlabel("TS")
            plt.ylabel("Probability Density")
            plt.yscale('log')
            plt.savefig(
                self.analysispath + '/TS_distribution.png',
                bbox_inches='tight')

    def __str__(self):
        r"""Print a message with info about the source once
        the analysis is running"""
        int_str = '*'*80
        int_str += '\n*' + ' '*78 + '*\n'
        int_str += '*' + ' '*((78-len(self.name))//2) + self.name + ' '*((78-len(self.name))//2 + len(self.name)%2) + '*'
        int_str += '\n*' + ' '*78 + '*\n'
        int_str += '*'*80 + '\n'
        int_str += '  '*5 + str(Time(self.start, format = 'mjd', scale = 'utc').iso)
        int_str += '  '*5 + str(Time(self.stop, format = 'mjd', scale = 'utc').iso) + '\n'
        return int_str

    def significance(self, p):
        r"""Given p value, report significance

        Parameters
        ------------
        p: float
            P value
            
        Returns
        --------
        sigma: float
            Significance
        """
        # sigma = stats.norm.isf(p)
        sigma = np.sqrt(2)*erfinv(1-2*p)
        return sigma

    def save_results(self, alt_path=None):
        r"""Once analysis has been performed, save
        all relevant info to a new directory.
        Saves to a file called [analysisid]_results.pickle

        Parameters
        -----------
        alt_path: string
            optional specified directory to save results to
        
        Returns
        -----------
        save_items: dict
            dictionary of saved analysis parameters        
        """
        if alt_path is None:
            print("Saving results to directory:\n\t{}".format(self.analysispath))
            with open(self.analysispath + '/' + self.analysisid + '_' + 'results.pickle', 'wb') as f:
                pickle.dump(self.save_items, f, protocol=pickle.HIGHEST_PROTOCOL)
            print("Results successfully saved")
            return self.save_items       
        else:
            print("Saving to specified directory")
            with open(alt_path + '/' + self.analysisid + '_' + 'results.pickle', 'wb') as f:
                pickle.dump(self.save_items, f, protocol=pickle.HIGHEST_PROTOCOL)
            print("Results successfully saved")
            return self.save_items

    def plot_ontime(self, with_contour=False, contour_files=None, label_events=False):
        r"""Plots ontime events on the full skymap and a 
        zoomed in version near the scan best-fit

        Parameters
        -----------
        with_contour: bool
            plots the 90% containment contour of a skymap (default False)
        contour_files: string
            text file containing skymap contours to be plotted (default None)
        label_events: bool
            adds a number label to events on skymap (default False)

        """
        
        try:
            self.plot_skymap_zoom(with_contour=with_contour, contour_files=contour_files)
        except Exception as e:
            print('Failed to make skymap zoom plot')

        try:
            self.plot_skymap(with_contour=with_contour, contour_files=contour_files, label_events=label_events) 
        except Exception as e:
            print('Failed to make FULL skymap plot')

    def plot_skymap_zoom(self, with_contour=False, contour_files=None):
        r"""Make a zoomed in portion of a skymap with
        all ontime neutrino events within a certain range
        Outputs a plot (in png and pdf formats) to the analysis path

        Parameters
        ------------
        with_contour: bool
            plots the 90% containment contour of a skymap (default False)
        contour_files: string
            text file containing skymap contours to be plotted (default None)

        """
        events = self.llh_exp # flattened array of exp per sample
        events = events[(events['time'] < self.stop) & (events['time'] > self.start)]

        col_num = 5000
        seq_palette = sns.color_palette("icefire", col_num)
        lscmap = mpl.colors.ListedColormap(seq_palette)

        rel_t = np.array((events['time'] - self.start) * col_num / (self.stop - self.start), dtype = int)
        cols = np.array([seq_palette[j] for j in rel_t])

        if self.skymap is not None:
            skymap = self.skymap
            ra = self.skymap_fit_ra
            dec = self.skymap_fit_dec
            label_str = "Scan Hot Spot"
            cmap = None
        else:
            skymap = np.zeros(hp.nside2npix(self._nside))
            ra = self.ra
            dec = self.dec
            label_str = self.name
            cmap = mpl.colors.ListedColormap([(1.,1.,1.)] * 50)

        plotting_utils.plot_zoom(skymap, ra, dec, "", range = [0,10], reso=3., cmap = cmap)

        if self.skipped is not None:
            try:
                msk = events['run'] == int(self.skipped[0][0])
                msk *= events['event'] == int(self.skipped[0][1])
                plotting_utils.plot_events(self.skipped_event['dec'], self.skipped_event['ra'], 
                    self.skipped_event['sigma']*self._angScale, 
                    ra, dec, 2*6, sigma_scale=1.0, constant_sigma=False, 
                    same_marker=True, energy_size=True, col = 'grey', 
                    with_dash=True)
                events = events[~msk]
                cols = cols[~msk]
            except:
                print("Removed event not in GFU")

        if (self.stop - self.start) <= 21.:
            plotting_utils.plot_events(events['dec'], events['ra'], events['sigma']*self._angScale, ra, dec, 2*6, sigma_scale=1.0,
                    constant_sigma=False, same_marker=True, energy_size=True, col = cols)
        else:
            #Long time windows means don't plot contours
            plotting_utils.plot_events(events['dec'], events['ra'], events['sigma']*self._angScale, ra, dec, 2*6, sigma_scale=None,
                    constant_sigma=False, same_marker=True, energy_size=True, col = cols)

        if contour_files is not None:
            cont_ls = ['solid', 'dashed']
            contour_counter = 0
            for c_file in contour_files:
                cont = np.loadtxt(c_file, skiprows=1)
                cont_ra = cont.T[0]
                cont_dec = cont.T[1]
                label = 'Millipede 50\%, 90\% (160427A syst.)' \
                    if contour_counter == 0 else ''
                hp.projplot(np.pi/2. - cont_dec, cont_ra, linewidth=3., 
                    color='k', linestyle=cont_ls[contour_counter], coord='C', 
                    label=label)
                contour_counter += 1
        
        if with_contour:
            probs = hp.pixelfunc.ud_grade(self.skymap, 64, power=-2)
            probs = probs/np.sum(probs)
            ### plot 90% containment contour of PDF
            levels = [0.9]
            theta, phi = plotting_utils.plot_contours(levels, probs)
            hp.projplot(theta[0], phi[0], linewidth=2., c='k')
            for i in range(1, len(theta)):
                hp.projplot(theta[i], phi[i], linewidth=2., c='k', label=None)

        plt.scatter(0,0, marker='*', c = 'k', s = 130, label = label_str) 

        plt.legend(loc = 2, ncol=1, mode = 'expand', fontsize = 18.5, framealpha = 0.95)
        plotting_utils.plot_color_bar(range=[0,6], cmap=lscmap, col_label=r"IceCube Event Time",
                    offset=-50, labels = [r'-$\Delta T \Bigg/ 2$', r'+$\Delta T \Bigg/ 2$'])
        plt.savefig(self.analysispath + '/' + self.analysisid + 'unblinded_skymap_zoom.png',bbox_inches='tight')
        plt.savefig(self.analysispath + '/' + self.analysisid + 'unblinded_skymap_zoom.pdf',bbox_inches='tight', dpi=300)
        plt.close()

    def plot_skymap(self, with_contour=False, contour_files=None, label_events=False):
        r""" Make skymap with event localization and all
        neutrino events on the sky within the given time window
        Outputs a plot in png format to the analysis path

        Parameters
        -----------
        with_contour: bool
            plots the 90% containment contour of a skymap (default False)
        contour_files: string
            text file containing skymap contours to be plotted (default None)
        label_events: bool
            adds a number label to events on skymap (default False)

        """
        events = self.llh_exp
        events = events[(events['time'] < self.stop) & (events['time'] > self.start)]

        col_num = 5000
        seq_palette = sns.color_palette("icefire", col_num)
        lscmap = mpl.colors.ListedColormap(seq_palette)

        rel_t = np.array((events['time'] - self.start) * col_num / (self.stop - self.start), dtype = int)
        cols = [seq_palette[j] for j in rel_t]

        # Set color map and plot skymap
        pdf_palette = sns.color_palette("Blues", 500)
        cmap = mpl.colors.ListedColormap(pdf_palette)
        cmap.set_under("w")

        if self.skymap is None:
            skymap = np.zeros(hp.nside2npix(self._nside))
            max_val = 1.
            cmap = mpl.colors.ListedColormap([(1.,1.,1.)] * 50)
        else:
            skymap = self.skymap
            max_val = max(skymap)
            min_val = min(skymap)
        moll_cbar = True if self.skymap is not None else None
        
        try:
            hp.mollview(skymap, coord='C', cmap=cmap, rot=180, cbar=moll_cbar) #cbar=None) #
        except Exception as e:
            if min_val<1.0e-16:
                #for some reason, sometimes have an underflow issue here
                skymap[skymap<1.0e-16] = 0.
            else: 
                print(e)
                print('Failed to make all-sky unblinded skymap. Retry making full skymap! ')
                return

        plt.text(2.0,0., r"$0^\circ$", ha="left", va="center")
        plt.text(1.9,0.45, r"$30^\circ$", ha="left", va="center")
        plt.text(1.4,0.8, r"$60^\circ$", ha="left", va="center")
        plt.text(1.9,-0.45, r"$-30^\circ$", ha="left", va="center")
        plt.text(1.4,-0.8, r"$-60^\circ$", ha="left", va="center")
        plt.text(2.0, -0.15, r"$0^\circ$", ha="center", va="center")
        plt.text(1.333, -0.15, r"$60^\circ$", ha="center", va="center")
        plt.text(.666, -0.15, r"$120^\circ$", ha="center", va="center")
        plt.text(0.0, -0.15, r"$180^\circ$", ha="center", va="center")
        plt.text(-.666, -0.15, r"$240^\circ$", ha="center", va="center")
        plt.text(-1.333, -0.15, r"$300^\circ$", ha="center", va="center")
        plt.text(-2.0, -0.15, r"$360^\circ$", ha="center", va="center")

        hp.graticule(verbose=False)

        theta=np.pi/2 - events['dec']
        phi = events['ra']

        #plot 90% containment errors
        sigma_90 = events['sigma']*self._angScale

        # plot events on sky with error contours
        handles=[]
        hp.projscatter(theta,phi,c=cols,marker='x',label='Event',coord='C', zorder=5)
        if label_events:
            for j in range(len(theta)):
                hp.projtext(theta[j], phi[j]-0.11, '{}'.format(j+1), color='red', fontsize=18, zorder=6)
        handles.append(Line2D([0], [0], marker='x', ls='None', label='Event'))

        if (self.stop - self.start) <= 0.5:        #Only plot contours if less than 2 days
            for i in range(events['ra'].size):
                my_contour = plotting_utils.contour(events['ra'][i], 
                                    events['dec'][i],sigma_90[i], self._nside)
                hp.projplot(my_contour[0], my_contour[1], linewidth=2., 
                                    color=cols[i], linestyle="solid",coord='C', zorder=5)

        if self.skymap is None:
            src_theta = np.pi/2. - self.dec
            src_phi = self.ra
            hp.projscatter(src_theta, src_phi, c = 'k', marker = '*',
                                label = self.name, coord='C', s=350)
            handles.append(Line2D([0], [0], marker='*', c='k', ls='None', label=self.name))

        if contour_files is not None:
            cont_ls = ['solid', 'dashed']
            contour_counter = 0
            for c_file in contour_files:
                cont = np.loadtxt(c_file, skiprows=1)
                cont_ra = cont.T[0]
                cont_dec = cont.T[1]
                hp.projplot(np.pi/2. - cont_dec, cont_ra, linewidth=3., 
                    color='k', linestyle=cont_ls[contour_counter], coord='C')
                contour_counter += 1

        if with_contour and self.skymap is not None:
            probs = hp.pixelfunc.ud_grade(self.skymap, 64, power=-2)
            probs = probs/np.sum(probs)
            ### plot 90% containment contour of PDF
            levels = [0.9]
            theta, phi = plotting_utils.plot_contours(levels, probs)
            hp.projplot(theta[0], phi[0], linewidth=2., c='k', label='Skymap (90\% cont.)')
            handles.append(Line2D([0], [0], lw=2, c='k', label=r"Skymap (90\% cont.)"))
            for i in range(1, len(theta)):
                hp.projplot(theta[i], phi[i], linewidth=2., c='k')
        
        plt.title(self.name.replace('_', ' '))
        plt.legend(loc=1, handles=handles)
        try: 
            plt.savefig(self.analysispath + '/' + self.analysisid + 'unblinded_skymap.png',bbox_inches='tight')
        except:
            plt.title('Fast Response Skymap')
            plt.savefig(self.analysispath + '/' + self.analysisid + 'unblinded_skymap.png',bbox_inches='tight')
        plt.close()

    def generate_report(self):
        r"""Generates report using class attributes
        and the ReportGenerator Class
        Makes full report and outputs as pdf
        """
        report = FastResponseReport(self)
        report.generate_report()
        report.make_pdf()
    

        
class PointSourceFollowup(FastResponseAnalysis):
    """
    Class for point-source or extended source followup
    i.e. there is a fixed location on the sky, not a healpy skymap

    The arguments index and dataset override the defaults of the class.
    """
    _nside = 256
    def __init__(self, name, ra, dec, tstart, tstop, extension=None,
                 index=None, dataset=None, fix_index=None,
                 skipped=None, outdir=None, save=True, seed=None):

        if index is not None:
            self._index = float(index)
        if dataset is not None:
            self._dataset = dataset
        if fix_index is not None:
            self._fix_index = fix_index
            self._float_index = not self._fix_index
        
        super().__init__(name, tstart, tstop, skipped=skipped, seed=seed,
                       outdir=outdir, save=save, extension=extension)
        

        

        self.ra = np.deg2rad(ra)
        self.dec = np.deg2rad(dec)
        
        if (self.ra < 0 or self.ra > 2*np.pi) or \
           (self.dec < -np.pi/2. or self.dec > np.pi/2.):
            print("Right ascension and declination are not valid")
            sys.exit()
        else:
            self.save_items['ra'] = self.ra
            self.save_items['dec'] = self.dec
        
        self.skymap, self.nside, self.ipix_90 = None, None, None

    def __str__(self):
        int_str = super().__str__()
        int_str += ' '*10 + 'RA, dec.: {:.2f}, {:+.3f}'.format(
            self.ra*180. / np.pi, self.dec*180. / np.pi)
        exts = 0. if self.extension is None else self.extension
        int_str += ' '*10 + 'Extension: {:.2f}'.format(exts * 180. / np.pi)
        int_str += '\n\n'
        return int_str

    def run_background_trials(self, ntrials=1000):
        r""" Method to run background trials for point source analysis

        Parameters
        -----------
        ntrials: int
            number of background trials to run (default 1000)
        """
        trials = self.llh.do_trials(
            int(ntrials),
            src_ra=self.ra,
            src_dec=self.dec)
        self.tsd = trials['TS']
        self.save_items['tsd'] = self.tsd

    def initialize_injector(self, e_range=(0., np.inf)):
        r"""Method to make relevant injector in Skylab.
        Used in calculating upper limit

        Parameters
        -----------
        e_range: tuple
            optional energy range allowed for injector. default (0., np.inf)

        Returns
        --------
        inj: Skylab injector object
            Point source injector using skylab
        """
        inj = PointSourceInjector(
            gamma = self.index, 
            E0 = 1000., 
            e_range=e_range)
        if self.multi:
            temporal_model = {enum:_llh.temporal_model for enum,_llh in self.llh._samples.items()}
        else:
            temporal_model = self.llh.temporal_model
        inj.fill(
            self.dec,
            self.llh.exp,
            self.llh.mc,
            self.llh.livetime,
            temporal_model=temporal_model)
        self.inj = inj
        self.save_items['E0'] = self.inj.E0

    def unblind_TS(self):
        r""" Unblind TS at one location for a point source

        Returns
        --------
        ts: float
            unblinded test statistic
        ns: float
            best fit ns
        """ 
        # Fix the case of getting best-fit gamma
        # TODO: What if gamma is floated
        ts, ns = self.llh.fit_source(src_ra=self.ra, src_dec=self.dec)
        params = ns.copy()
        params.pop('nsignal')
        self.ns_params = params
        ns = ns['nsignal']
        if self._verbose:
            print("TS = {}".format(ts))
            print("ns = {}".format(ns))
            for par, val in params.items():
                print(f"{par} = {val:.3f}")
            print("\n\n")
        self.ts, self.ns = ts, ns
        self.save_items['ts'] = ts
        self.save_items['ns'] = ns
        # TODO alternatively change ReportGenerator to report ns_params generically
        if 'gamma' in params:
            self.gamma = params['gamma']
            for par, val in params.items():
                if par in self.save_items:
                    if self._verbose:
                        print(f'Warning, not saving {par} as save_items already has such a key')
                self.save_items.setdefault(par, val)
        return ts, ns

    def find_coincident_events_single(self, llh):
        r"""Find "coincident events" for the analysis.
        These are ontime events that have a spatial times energy weight greater than 10
        """
        # FIXME these will not work for Multi
        # either wrap class or modify attributes

        spatial_weights = llh.llh_model.signal(
            self.ra, self.dec, llh._events, 
            src_extension=self.extension)[0] / llh._events['B']
        energy_ratio, _ = llh.llh_model.weight(
           llh._events, **self.ns_params)
        temporal_weights =llh.temporal_model.signal(llh._events)
        msk = spatial_weights * energy_ratio * temporal_weights > 10
        coincident_events = []
        if len(spatial_weights[msk]) > 0:
            coincident_events = []
            for ev, s_w, en_w in zip(llh._events[msk], 
                        spatial_weights[msk], energy_ratio[msk]):
                coincident_events.append({})
                for key in ['run', 'event', 'ra', 'dec', 'sigma', 'logE', 'time']:
                    coincident_events[-1][key] = ev[key]
                del_psi = sensitivity_utils.deltaPsi([ev['dec']], [ev['ra']], [self.dec], [self.ra])[0]
                coincident_events[-1]['delta_psi'] = del_psi 
                coincident_events[-1]['spatial_w'] = s_w
                coincident_events[-1]['energy_w'] = en_w
        return coincident_events
    
    def find_coincident_events(self):
        r"""Find "coincident events" for the analysis.
        These are ontime events that have a spatial times energy weight greater than 10.
        These will be combined from all used samples, with a key 'enum' to distinguish.
        """
        # TODO clarify the docstring: actually it includes the temporal weight!
        coincident_events = []
        if self.multi:
            for enum, _sam in self.llh._samples.items():
                _coincident_events = self.find_coincident_events_single(_sam)
                for _ev in _coincident_events:
                    _ev['enum'] = enum
                coincident_events += _coincident_events
                
        else:
            coincident_events = self.find_coincident_events_single(self.llh)
        self.coincident_events = coincident_events
        self.save_items['coincident_events'] = coincident_events

    def ns_scan(self, params = None): 
        r""" Calculate the llh as a function of ns, to 
        see how shallow the likelihood space is.

        Parameters
        -----------
        params: mappable
            params dictionary for skylab.ps_llh object. Default is 
            that of this analysis (set below, using self.spectrum)
            
        Returns
        ---------
        xs, delta_llh: array-like
            Values of ns scanned and corresponding -2*delta_llh values"""

        if params is None:
            if self._fix_index:
                params = {'spectrum': str(self.spectrum)}
            else:
                params = self.ns_params
                # TODO or minimize per point?

        bounds = self.ns * 2.5
        if bounds == 0.0:
            bounds = 1.0
        xs = np.linspace(0., bounds, 101)
        llhs = []
        for n in xs:
            llhs.append(self.llh.llh(n, **params)[0])
        llhs = np.array(llhs)
        best_llh = self.llh.llh(self.ns, **params)[0]
        delta_llh = -2.*(llhs - best_llh)

        if self.save_output:
            fig, ax = plt.subplots()
            plt.plot(xs, delta_llh, lw = 3)
            plt.axvline(self.ns, ls = '--', color = 'gray', lw = 1)
            plt.xlabel(r'$n_s$', fontsize = 18)
            plt.ylabel(r'$-2\times \Delta \mathrm{LLH}$', fontsize = 18)
            plt.savefig(self.analysispath + '/llh_ns_scan.png', bbox_inches='tight', dpi=200)

        self.ns_profile = (xs, delta_llh)
        self.save_items['ns_scan'] = xs, delta_llh
        return xs, delta_llh

    def upper_limit(self, n_per_sig=100, p0=None):
        r"""After calculating TS, find upper limit
        Assuming an E^-2 spectrum

        Parameters
        -----------
        n_per_sig: int
            number of trials per injected signal events to run
        p0: None, scalar, or n-length sequence
            Initial guess for the parameters passed to curve_fit (optional, default None)

        Returns
        --------
        upperlimit: float
            Value of E^2 dN / dE in units of TeV / cm^2 
        """
        # TODO this docstring is not general enough,
        # and actually the return value is from mu2flux but does not match its units
        # because injector index is set by self._index, and for multisample we want softer indices.
        if self.inj is None:
            self.initialize_injector()
        ninj = np.array([1., 1.5, 2., 2.5, 3., 4., 5., 6.])
        passing = []
        for n in ninj:
            results = self.llh.do_trials(
                n_per_sig, src_ra=self.ra, src_dec=self.dec,
                injector=self.inj, mean_signal=n, poisson=True)
            msk = results['TS'] > self.ts
            npass = len(results['TS'][msk])
            passing.append((n, npass, n_per_sig))
            
        signal_fluxes, passing, number = list(zip(*passing))
        signal_fluxes = np.array(signal_fluxes)
        passing = np.array(passing, dtype=float)
        number = np.array(number, dtype=float)
        passing /= number
        errs = sensitivity_utils.binomial_error(passing, number)
        fits, plist = [], []
        for func, func_name in [(sensitivity_utils.chi2cdf, 'chi2'),
                                (sensitivity_utils.erfunc, 'Error function'),
                                (sensitivity_utils.incomplete_gamma, 
                                'Incomplete gamma')]:
            try:
                fits.append(sensitivity_utils.sensitivity_fit(
                    signal_fluxes, passing, errs, func, p0=p0))
                plist.append(fits[-1]['pval'])
            except:
                print(f"{func_name} fit failed in upper limit calculation")
        #Find best fit of the three, make it look different in plot
        plist = np.array(plist)
        best_fit_ind = np.argmax(plist)
        fits[best_fit_ind]['ls'] = '-'
        self.upperlimit = self.inj.mu2flux(fits[best_fit_ind]['sens'])
        self.upperlimit_ninj = fits[best_fit_ind]['sens']
        # E^2 dN/dE at self.inj.E0
        upperlimit_fluence = self.upperlimit * self.duration * 86400. * self.inj.E0**2

        fig, ax = plt.subplots()
        for fit_dict in fits:
            lw = 1. if fit_dict['ls'] == '-' else 1.5
            ax.plot(
                fit_dict['xfit'], fit_dict['yfit'], 
                label = r'{}: $\chi^2$ = {:.2f}, d.o.f. = {}'.format(
                    fit_dict['name'], fit_dict['chi2'], fit_dict['dof']),
                ls = fit_dict['ls'],
                lw=lw)
            if fit_dict['ls'] == '-':
                ax.axhline(0.9, color = 'm', linewidth = 0.3, linestyle = '-.')
                ax.axvline(fit_dict['sens'], color = 'm', linewidth = 0.3, linestyle = '-.')
                limit_annotation = r'Sens. = {:.2f} events'.format(fit_dict['sens']) + '\n'
                limit_annotation +=  r' = {:.1e}'.format(upperlimit_fluence) + r' GeV cm$^{-2}$' + '\n'
                #ax.text(3.5, 0.8, r'Sens. = {:.2f} events'.format(fit_dict['sens']), fontsize = 16)
                #ax.text(3.5, 0.7, r' = {:.1e}'.format(upperlimit_fluence) + r' GeV cm$^{-2}$', fontsize = 16)
                if self.index != 2:
                    # E^2 F not constant in energy, state pivot energy in label
                    #ax.text(3.5, 0.7, r' = {:.1e}'.format(upperlimit_fluence) + r' GeV cm$^{-2}$' + f'\nat {self.inj.E0:.0f} GeV', fontsize = 16)
                    limit_annotation += f'at {self.inj.E0:.0f} GeV'
                #ax.text(5, 0.5, r'Sens. = {:.2e}'.format(self.inj.mu2flux(fit_dict['sens'])) + ' GeV^-1 cm^-2 s^-1')
                ax.annotate(limit_annotation, (0.6, 0.8), ha = 'left', va = 'top', xycoords = 'axes fraction', fontsize = 16)
        ax.errorbar(signal_fluxes, passing, yerr=errs, capsize = 3, linestyle='', marker = 's', markersize = 2)
        ax.legend(loc=4, fontsize = 14)
        ax.set_xlabel(r'$\langle n_{inj} \rangle$', fontsize = 14)
        ax.set_ylabel(r'Fraction TS $>$ unblinded TS', fontsize = 14)
        if self.save_output:
            plt.savefig(self.analysispath + '/upper_limit_distribution.png', bbox_inches='tight', dpi=200)
        print("Found upper limit of {:.2f} events".format(self.upperlimit_ninj))
        self.save_items['upper_limit'] = self.upperlimit
        self.save_items['upper_limit_ninj'] = self.upperlimit_ninj
        return self.upperlimit

    def make_dNdE_single(self):
        r"""Make an E^-2 or E^-2.5 dNdE with the central 90% 
        for the most relevant declination band 
        (+/- 5 deg around source dec)
        """
        dec_mask_1 = self.llh.mc['dec'] > self.dec - (5. * np.pi / 180.)
        dec_mask_2 = self.llh.mc['dec'] < self.dec + (5. * np.pi / 180.)
        dec_mask_3, dec_mask_4 = None, None
        dec_mask = dec_mask_1 * dec_mask_2
        fig, ax = plt.subplots(figsize = (8,5))
        fig.set_facecolor('white')
        lab = None
        delta_gamma = -1. * self.index + 1.
        a = plt.hist(self.llh.mc['trueE'][dec_mask], bins = np.logspace(1., 8., 50), 
                weights = self.llh.mc['ow'][dec_mask] * np.power(self.llh.mc['trueE'][dec_mask], delta_gamma) / self.llh.mc['trueE'][dec_mask], 
                histtype = 'step', linewidth = 2., color = sns.xkcd_rgb['windows blue'], label = lab)
        plt.yscale('log')
        plt.xscale('log')
        plt.grid(which = 'major', alpha = 0.25)
        plt.xlabel('Energy (GeV)', fontsize = 24)
        cdf = np.cumsum(a[0]) / np.sum(a[0])
        low_5 = np.interp(0.05, cdf, a[1][:-1])
        median = np.interp(0.5, cdf, a[1][:-1])
        high_5 = np.interp(0.95, cdf, a[1][:-1])
        self.low5 = low_5
        self.high5 = high_5
        plt.axvspan(low_5, high_5, color = sns.xkcd_rgb['windows blue'], alpha = 0.25, label="Central 90\%")
        lab = 'Median'
        plt.axvline(median, c = sns.xkcd_rgb['windows blue'], alpha = 0.75, label = lab)
        plt.xlim(min(min(self.low5), 1e1), max(max(self.high5), 1e8))
        plt.legend(loc=4, fontsize=18)
        plt.savefig(self.analysispath + '/central_90_dNdE.png',bbox_inches='tight')

        print('90% Central Energy Range: {}, {} GeV'.format(round(low_5), round(high_5)))
        self.save_items['energy_range'] = (self.low5, self.high5)

    def make_dNdE_multi(self):
        r"""Make an E^-2 or E^-2.5 dNdE with the central 90% 
        for the most relevant declination band 
        (+/- 5 deg around source dec)
        """
        low5 = []
        high5 = []
        fig, ax = plt.subplots(figsize = (8,5))
        fig.set_facecolor('white')
        for enum in self.llh._samples:
            llh = self.llh._samples[enum]
            dataset = self.dataset[enum].replace('_', ' ')

            dec_mask_1 = llh.mc['dec'] > self.dec - (5. * np.pi / 180.)
            dec_mask_2 = llh.mc['dec'] < self.dec + (5. * np.pi / 180.)
            dec_mask_3, dec_mask_4 = None, None
            dec_mask = dec_mask_1 * dec_mask_2
            
            
            lab = dataset
            color = f'C{enum+1}'
            delta_gamma = -1. * self.index + 1.
            a = plt.hist(llh.mc['trueE'][dec_mask], bins = np.logspace(1., 8., 50), 
                    weights = llh.mc['ow'][dec_mask] * np.power(llh.mc['trueE'][dec_mask], delta_gamma) / llh.mc['trueE'][dec_mask], 
                    histtype = 'step', linewidth = 2., color = color, label = lab)
            cdf = np.cumsum(a[0]) / np.sum(a[0])
            low_5 = np.interp(0.05, cdf, a[1][:-1])
            median = np.interp(0.5, cdf, a[1][:-1])
            high_5 = np.interp(0.95, cdf, a[1][:-1])
        
            plt.axvspan(low_5, high_5, color = color, alpha = 0.25, label="Central 90\%")
            lab = 'Median'
            plt.axvline(median, c = color, alpha = 0.75, label = lab)
            low5.append(low_5)
            high5.append(high_5)
            print('{} 90% Central Energy Range: {}, {} GeV'.format(dataset, round(low_5), round(high_5)))
            
        plt.yscale('log')
        plt.xscale('log')
        plt.grid(which = 'major', alpha = 0.25)
        plt.xlabel('Energy (GeV)', fontsize = 24)

        plt.xlim(1e1, 1e8)
        plt.legend(loc='upper left',  bbox_to_anchor=(1,1), fontsize=18)
        plt.savefig(self.analysispath + '/central_90_dNdE.png',bbox_inches='tight')

        self.low5 = low5
        self.high5 = high5
        self.save_items['energy_range'] = tuple(zip(self.low5, self.high5))

    def make_dNdE(self):
        if self.multi:
            self.make_dNdE_multi()
        else:
            self.make_dNdE_single()

    def write_circular(self):
        pass

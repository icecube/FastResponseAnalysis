r''' General Fast Response Analysis Class.

    Author: Alex Pizzuto
    Date: 2021
    '''

from abc import abstractmethod
import os, sys, time, subprocess
import pickle, dateutil.parser, logging, warnings

import h5py
import healpy               as hp
import numpy                as np
import seaborn              as sns
import matplotlib           as mpl
import matplotlib.pyplot    as plt
from astropy.time           import Time
from scipy.special          import erfinv
from matplotlib.lines       import Line2D

from skylab.datasets        import Datasets
from skylab.llh_models      import EnergyLLH
from skylab.priors          import SpatialPrior
from skylab.ps_injector     import PointSourceInjector
from skylab.ps_llh          import PointSourceLLH
from skylab.ps_injector     import PriorInjector
from skylab.spectral_models import PowerLaw 
from skylab.temporal_models import BoxProfile, TemporalModel
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
    _season_names = [f"IC86, 201{y}" for y in range(1, 10)]
    _nb_days = 10
    _ncpu = 5

    def __init__(self, name, tstart, tstop,
                 skipped=None, seed=None,
                 outdir=None, save=True,
                 extension=None,
                 index=None,
                 fix_index=None,
                 dataset=None,
                 ):
        self.name = name

        if index is not None:
            self._index = float(index)
        if dataset is not None:
            self._dataset = dataset
        if fix_index is not None:
            self._fix_index = fix_index
            self._float_index = not self._fix_index
        
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
        self.llh = self.initialize_llh(skipped=skipped, scramble=self.scramble)
        self.inj = None

    @property
    def dataset(self):
        """Returns the dataset used"""
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
        """Returns a flat array of experimental data loaded into the LLH,
        limited to fields used for plotting and amended with an enum=0."""
        if not hasattr(self, '_llh_exp'):
            exp = self.llh.exp
            merged_dtype = np.dtype([('run', '<i4'), ('event', '<i4'), ('time', '<f4'), ('ra', '<f4'), ('dec', '<f4'), ('sigma', '<f4'), ('enum', '<i4')])
            exp = rf.drop_fields(exp, [_field for _field in exp.dtype.names if _field not in merged_dtype.names])
            enum = 0
            exp = rf.append_fields(exp, 'enum', np.full( exp.size, enum))
            exp = exp.astype(merged_dtype)
            self._llh_exp = exp
        return self._llh_exp
    

    def get_data(self, livestream_start=None, livestream_stop=None):
        """
        Gets the skylab data and MC from querying the i3live livestream

        Parameters
        -----------
        livestream_start: float
            (optional) start time for getting data (MJD) from i3Live. Default is start time of the analysis - 5 days
        livestream_stop: float
            (optional) stop time for getting data (MJD). Default is stop time of the analysis
        """
        if self._verbose:
            print("Grabbing data")

        dset = Datasets[self.dataset]
        #if self.stop < 58933.0: 
        if self.stop < 59215:
            if self._verbose:
                print("Old times, just grabbing archival data")
            exps, grls = [], []
            for season in self._season_names:
                exp, mc, livetime = dset.season(season, floor=self._floor)
                grl = dset.grl(season)
                exps.append(exp)
                grls.append(grl)
            if self.stop > 58933.0:
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
        else:
            if self._verbose:
                print("Recent time: querying the i3live database")
            if livestream_start is None or livestream_stop is None:
                livestream_start = self.start - 6.
                livestream_stop = self.stop
            exp, mc, livetime, grl = dset.livestream(
                livestream_start, livestream_stop,
                append=self._season_names, 
                floor=self._floor)
        exp.sort(order='time')
        grl.sort(order='run')
        livetime = grl['livetime'].sum()

        # workaround while not all datasets have livestream yet
        if 'livestream' in dset.season_names():
            reference_season = 'livestream'
        else:
            reference_season = self._season_names[0]
        sinDec_bins = dset.sinDec_bins(reference_season)
        energy_bins = dset.energy_bins(reference_season)

        self.exp = exp
        self.mc = mc
        self.grl = grl
        self.livetime = livetime
        self.dset = dset
        self.sinDec_bins = sinDec_bins
        self.energy_bins = energy_bins
        
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

        if self._fix_index:
            llh_model = EnergyLLH(
                twodim_bins=[self.energy_bins, self.sinDec_bins],
                allow_empty=True,
                spectrum=PowerLaw(A=1, gamma=self.index, E0=1000.),
                ncpu=self._ncpu) 
        elif self._float_index:
            llh_model = EnergyLLH(
                twodim_bins=[self.energy_bins, self.sinDec_bins],
                allow_empty=True,
                bounds=self._index_range,
                seed = self.index,
                ncpu=self._ncpu)
        
        box = TemporalModel(
            grl=self.grl,
            poisson_llh=True,
            days=self._nb_days,
            signal=BoxProfile(self.start, self.stop))
        
        if skipped is not None:
            self.exp = self.remove_event(self.exp, self.dset, skipped)

        if self.extension is not None:
            src_extension = float(self.extension)
            self.save_items['extension'] = src_extension
        else:
            src_extension = None

        llh = PointSourceLLH(
            self.exp,                      # array with data events
            self.mc,                       # array with Monte Carlo events
            self.livetime,                 # total livetime of the data events
            ncpu=self._ncpu,               # use 10 CPUs when computing trials
            scramble=scramble,             # set to False for unblinding
            timescramble=True,             # not just RA scrambling
            llh_model=llh_model,           # likelihood model
            temporal_model=box,            # use box for temporal model
            nsource_bounds=(0., 1e3),      # bounds on fitted ns
            src_extension = src_extension, # Model symmetric extended source
            nsource=1.,                    # seed for nsignal fit
            seed=self.llh_seed)           

        return llh
    
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

    def plot_skymap_zoom(self, with_contour=False, contour_files=None, reso=3.):
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
        
        events = self.llh_exp
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

        plotting_utils.plot_zoom(skymap, ra, dec, "", range = [0,10], reso=reso, cmap = cmap)
        

        if self.skipped is not None:
            try:
                msk = events['run'] == int(self.skipped[0][0])
                msk *= events['event'] == int(self.skipped[0][1])
                # TODO here we don't know which sample the skipped event came from...
                plotting_utils.plot_events(self.skipped_event['dec'], self.skipped_event['ra'], 
                    self.skipped_event['sigma']*self._angScale, 
                    ra, dec, 2*6, sigma_scale=1.0, constant_sigma=False, 
                    same_marker=True, energy_size=True, col = 'grey', 
                    with_dash=True)
                events = events[~msk]
                cols = cols[~msk]
            except:
                print("Removed event not in GFU")
                # TODO print the actual dataset name
                # (once one implemention Multi code is merged)

        if (self.stop - self.start) <= 21.:
            sigma_scale = 1.0
        else:
            #Long time windows means don't plot contours
            sigma_scale = None

        for enum in np.unique(events['enum']):
            _mask = events['enum'] == enum
            _events = events[_mask]
            print(_events.size, self.datasets[enum], plotting_utils.skymap_style[enum])
            plotting_utils.plot_events(_events['dec'], _events['ra'], _events['sigma']*self._angScale,
                ra, dec, 2*6,
                sigma_scale=reso/3.,
                constant_sigma=False, same_marker=True, energy_size=True,
                col = cols[_mask], kw_style=plotting_utils.skymap_style[enum],
                )

        # plotting_utils.plot_events(events['dec'], events['ra'], events['sigma']*self._angScale,
        #         ra, dec, 2*6,
        #         sigma_scale=sigma_scale,
        #         constant_sigma=False, same_marker=True, energy_size=True, col = cols)

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
        #plt.close()

    def plot_skymap(self, events=None, with_contour=False, contour_files=None, label_events=False, label='GFU Event'):
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

        if events is None:
            events = self.llh.exp
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
        hp.projscatter(theta,phi,c=cols,marker='x',label=label,coord='C', zorder=5)
        if label_events:
            for j in range(len(theta)):
                hp.projtext(theta[j], phi[j]-0.11, '{}'.format(j+1), color='red', fontsize=18, zorder=6)
        handles.append(Line2D([0], [0], marker='x', ls='None', label=label))

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
    

class PriorFollowup(FastResponseAnalysis):
    """
    Class for skymap-based analyses
    """
    _pixel_scan_nsigma = 4.0
    _containment = 0.99
    _allow_neg = False
    _nside = 256

    def __init__(self, name, skymap_path, tstart, tstop, skipped=None, seed=None,
                 outdir=None, save=True, extension=None):

        super().__init__(name, tstart, tstop, skipped=skipped, seed=seed,
                       outdir=outdir, save=save, extension=extension)
        
        self.skymap_path = skymap_path
        try:
            skymap, skymap_header = hp.read_map(skymap_path, h=True, verbose=False)
            self.skymap_header = skymap_header
        except (TypeError, OSError):
            hdf_data = h5py.File(skymap_path, 'r')
            probs = hdf_data['PROBDENSITY'][()]
            area = 4*np.pi/probs.size
            probs *= area
            skymap = hp.pixelfunc.ud_grade(
                probs, self._nside, power=-2,
                order_in='NESTED', order_out='RING'
            )
            skymap_header = None
        skymap = self.format_skymap(skymap)
        self.skymap = skymap
        self.nside = hp.pixelfunc.get_nside(self.skymap)
        self.ipix_90 = self.ipixs_in_percentage(0.9)
        self.ra, self.dec, self.extension = None, None, extension
        self.save_items['skymap'] = skymap

    def __str__(self):
        int_str = super().__str__()
        int_str += ' '*10 + 'Skymap file:' + self.skymap_path
        int_str += '\n\n'
        return int_str

    def format_skymap(self, skymap):
        r"""Method to up or downgrade nside of a skymap to 
        the nside used in the analysis

        Parameters
        -----------
        skymap: array
            Healpix skymap from a file
        
        Returns
        ----------
        skymap: array
            Healpix skymap, with correct nside for use in FRA
        """
        if hp.pixelfunc.get_nside(skymap) != self._nside:
            skymap = hp.pixelfunc.ud_grade(skymap, self._nside, power=-2)
            skymap = skymap/skymap.sum()
        return skymap

    def initialize_injector(self, e_range=(0., np.inf)):
        r"""Method to make relevant injector in Skylab, done for analysis
        checks as well as for calculating sensitivities

        Parameters
        -----------
        e_range: tuple
            optional energy range allowed for injector. default (0., np.inf)

        Returns
        --------
        inj: Skylab injector object
            Spatial prior injector using skylab
        """

        print("Initializing Prior Injector")
        spatial_prior = SpatialPrior(self.skymap, containment = self._containment, allow_neg=self._allow_neg)
        self.spatial_prior = spatial_prior
        inj = PriorInjector(
            spatial_prior, 
            gamma=self.index, 
            e_range = e_range, 
            E0=1000., 
            seed = self.llh_seed)
        inj.fill(
            self.llh.exp,
            self.llh.mc,
            self.llh.livetime,
            temporal_model=self.llh.temporal_model)
        self.inj = inj

    def run_background_trials(self, ntrials=1000):
        r""" helper method to run background trials for the spatial prior analyses
        Note: this method is not used in running the analysis in realtime

        Returns
        --------
        ntrials: int
            number of trials to run (default 1000)
        """ 
        tsd = []
        spatial_prior = SpatialPrior(self.skymap, containment = self._containment, allow_neg=self._allow_neg)

        self.llh = self.initialize_llh(
            skipped = self.skipped, 
            scramble=True)

        toolbar_width = 50
        sys.stdout.write("[%s]" % (" " * toolbar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (toolbar_width+1))
        for i in range(int(ntrials)):
            if ntrials > 50 and i % (ntrials / 50) == 0:
                sys.stdout.write("#")
                sys.stdout.flush()
            val = self.llh.scan(
                0.0, 0.0, scramble = True, seed=i, 
                spatial_prior=spatial_prior,
                time_mask = [self.duration/2., self.centertime], 
                pixel_scan=[self.nside, self._pixel_scan_nsigma])
            try:
                tsd.append(val['TS_spatial_prior_0'].max())
            except ValueError:
                tsd.append(0.)
        tsd = np.asarray(tsd, dtype=float)
        self.tsd = tsd
        self.save_items['tsd'] = tsd

    def find_coincident_events(self):
        r"""Find coincident events for a skymap
        based analysis. These are ontime events that are also in the 
        90% contour of the skymap
        """
        t_mask=(self.llh.exp['time']<=self.stop)&(self.llh.exp['time']>=self.start)
        events = self.llh.exp[t_mask]
        exp_theta = 0.5*np.pi - events['dec']
        exp_phi   = events['ra']
        exp_pix   = hp.ang2pix(self.nside, exp_theta, exp_phi)
        overlap   = np.isin(exp_pix, self.ipix_90)
        events = events[overlap]

        # print nearby events, as a check (if needed)
        # msk1 = (self.llh.exp[t_mask]['ra'] < (self.skymap_fit_ra+np.radians(5)))*(self.llh.exp[t_mask]['ra'] > (self.skymap_fit_ra-np.radians(5)))
        # msk2 = (self.llh.exp[t_mask]['dec'] < (self.skymap_fit_dec+np.radians(5)))*((self.llh.exp[t_mask]['dec'] > self.skymap_fit_dec-np.radians(5)))
        # msk3 = msk1*msk2
        # print('Nearby events:')
        # print("[run, event, ra, dec, sigma, logE, time]")
        # for e in self.llh.exp[t_mask][msk3]: print([e[k] for k in ['run', 'event', 'ra', 'dec', 'sigma', 'logE', 'time']])

        if len(events) == 0:
            coincident_events = []
        else:
            coincident_events = []
            for ev in events:
                coincident_events.append({})
                for key in ['run', 'event', 'ra', 'dec', 'sigma', 'logE', 'time']:
                    coincident_events[-1][key] = ev[key]
                coincident_events[-1]['delta_psi'] = np.nan
                coincident_events[-1]['spatial_w'] = np.nan
                coincident_events[-1]['energy_w'] = np.nan
                coincident_events[-1]['in_contour'] = True
        self.coincident_events = coincident_events
        self.save_items['coincident_events'] = coincident_events

    def unblind_TS(self, custom_events=None):
        r""" Unblind TS, either sky scan for spatial prior,
        or just at one location for a point source

        Parameters
        -----------
        custom_events: array
            specific events for use in scan (UNUSED)
        
        Returns
        -----------
        ts: float
            unblinded test statistic
        ns: float
            unblinded fitted number of signal events at the best-fit location
        gamma: float
            unblinded fitted spectral index at the best-fit location
            returned only if spectral index is floated
        
        """ 
        t0 = time.time()
        spatial_prior = SpatialPrior(
            self.skymap, containment = self._containment,
            allow_neg=self._allow_neg)
        pixels = np.arange(len(self.skymap))
        t1 = time.time()
        print("Starting scan")
        val = self.llh.scan(
            0.0,0.0, scramble = False, spatial_prior=spatial_prior,
            time_mask = [self.duration/2., self.centertime],
            pixel_scan=[self.nside, self._pixel_scan_nsigma]
        )
        self.ts_scan = val
        t2 = time.time()
        print("finished scan, took {} s".format(t2-t1))
        try:
            ts = val['TS_spatial_prior_0'].max()
            max_prior = np.argmax(val['TS_spatial_prior_0'])
            ns = val['nsignal'][max_prior]
            if self._float_index:
                gamma_fit = val['gamma'][max_prior]
            else:
                gamma_fit = np.nan
            if ts != 0.0:
                self.skymap_fit_ra = val['ra'][max_prior]
                self.skymap_fit_dec = val['dec'][max_prior]
            else:
                max_pix = np.argmax(self.skymap)
                theta, phi = hp.pix2ang(self.nside, max_pix)
                dec = np.pi/2. - theta
                self.skymap_fit_ra = phi
                self.skymap_fit_dec = dec
            self.scanned_pixels = hp.ang2pix(
                self.nside, np.pi/2. - val['dec'], val['ra']
            )
        except Exception as e:
            print(e)
            ts, ns = 0., 0.
            if self._float_index:
                gamma_fit = 2.0
            else:
                gamma_fit = np.nan
            max_pix = np.argmax(self.skymap)
            theta, phi = hp.pix2ang(self.nside, max_pix)
            dec = np.pi/2. - theta
            self.skymap_fit_ra = phi
            self.skymap_fit_dec = dec
            self.scanned_pixels = None
        print("TS = {}".format(ts))
        print("ns = {}".format(ns))
        self.ts, self.ns = ts, ns
        if self._float_index:
            self.gamma = gamma_fit
            self.save_items['gamma'] = gamma_fit
            print(f'gamma = {gamma_fit}')
        print('\n\n')
        self.save_items['ts'] = ts
        self.save_items['ns'] = ns
        self.save_items['fit_ra'] = self.skymap_fit_ra
        self.save_items['fit_dec'] = self.skymap_fit_dec
        if self._float_index:
            return ts, ns, gamma_fit
        else:
            return ts, ns

    def upper_limit(self):
        """ UPPER LIMIT WITH SPATIAL PRIOR NOT YET IMPLEMENTED
        """
        print("Upper limit with spatial prior not yet implemented")
        pass

    def ipixs_in_percentage(self, percentage):
        """Finding ipix indices confined in a given percentage.
        
        Parameters
        ----------------
        skymap: ndarray
            array of probabilities for a skymap
        percentage : float
            fractional percentage from 0 to 1  

        Returns
        ---------
        ipix : numpy int array
            indices of pixels within percentage containment
        """
        # TODO: Find a more efficient way to do this
        skymap = self.skymap
        sort = sorted(skymap, reverse = True)
        cumsum = np.cumsum(sort)
        index, value = min(enumerate(cumsum),key=lambda x:abs(x[1]-percentage))

        index_hpx = list(range(0,len(skymap)))
        hpx_index = np.c_[skymap, index_hpx]

        sort_2array = sorted(hpx_index,key=lambda x: x[0],reverse=True)
        value_contour = sort_2array[0:index]

        j = 1 
        table_ipix_contour = [ ]
        for i in range (0, len(value_contour)):
            ipix_contour = int(value_contour[i][j])
            table_ipix_contour.append(ipix_contour)
    
        ipix = table_ipix_contour
          
        return np.asarray(ipix,dtype=int)

    def dec_skymap_range(self):
        r""" Compute minimum and maximum declinations within
        of the 90% contour of a given skymap

        Returns
        --------
        low: float
            lowest declination in 90%
        high: float
            highest declination in 90%
        """

        src_theta, src_phi = hp.pix2ang(self.nside, self.ipix_90)
        src_dec = np.pi/2. - src_theta
        src_dec = np.unique(src_dec)

        low = src_dec.min()
        high = src_dec.max()

        return low, high

    def make_dNdE(self):
        r"""Make an E^-2 or E^-2.5 dNdE with the central 90% 
        for the minimum and maximum declinations on the skymap
        """

        min_dec, max_dec = self.dec_skymap_range()
        dec_mask_1 = self.llh.mc['dec'] > min_dec - (5. * np.pi / 180.)
        dec_mask_2 = self.llh.mc['dec'] < min_dec + (5. * np.pi / 180.)
        dec_mask_3 = self.llh.mc['dec'] > max_dec - (5. * np.pi / 180.)
        dec_mask_4 = self.llh.mc['dec'] < max_dec + (5. * np.pi / 180.)
        dec_mask = dec_mask_1 * dec_mask_2

        # Start with the lower declination
        fig, ax = plt.subplots(figsize = (8,5))
        fig.set_facecolor('white')
        lab = 'Min dec.'
        delta_gamma = -1. * self.index + 1.
        a = plt.hist(self.llh.mc['trueE'][dec_mask], bins = np.logspace(1., 8., 50), 
                weights = self.llh.mc['ow'][dec_mask] * np.power(self.llh.mc['trueE'][dec_mask], delta_gamma) / self.llh.mc['trueE'][dec_mask], 
                histtype = 'step', linewidth = 2., color = sns.xkcd_rgb['windows blue'], label = lab)
        plt.yscale('log')
        plt.xscale('log')
        plt.grid(which = 'major', alpha = 0.25)
        plt.xlabel('Energy (GeV)', fontsize = 24)
        cdf = np.cumsum(a[0]) / np.sum(a[0])
        low_5_min_dec = np.interp(0.05, cdf, a[1][:-1])
        median_min_dec = np.interp(0.5, cdf, a[1][:-1])
        high_5_min_dec = np.interp(0.95, cdf, a[1][:-1])
        plt.axvspan(low_5_min_dec, high_5_min_dec, color = sns.xkcd_rgb['windows blue'], alpha = 0.25, label="Central 90\%")
        lab = 'Median (min dec.)'
        plt.axvline(median_min_dec, c = sns.xkcd_rgb['windows blue'], alpha = 0.75, label = lab)

        # Now do the same for the higher declinations
        dec_mask = dec_mask_3 * dec_mask_4
        a = plt.hist(self.llh.mc['trueE'][dec_mask], bins = np.logspace(1., 8., 50), 
            weights = self.llh.mc['ow'][dec_mask] * np.power(self.llh.mc['trueE'][dec_mask], delta_gamma) / self.llh.mc['trueE'][dec_mask], 
            histtype = 'step', linewidth = 2., color = sns.xkcd_rgb['dark navy blue'], label = 'Max dec.')
        cdf = np.cumsum(a[0]) / np.sum(a[0])
        low_5_max_dec = np.interp(0.05, cdf, a[1][:-1])
        median_max_dec = np.interp(0.5, cdf, a[1][:-1])
        high_5_max_dec = np.interp(0.95, cdf, a[1][:-1])
        plt.axvspan(low_5_max_dec, high_5_max_dec, color = sns.xkcd_rgb['dark navy blue'], alpha = 0.25)
        plt.axvline(median_max_dec, c = sns.xkcd_rgb['dark navy blue'], alpha = 0.75, label = "Median (max dec.)", ls = '--')
        plt.xlim(1e1, 1e8)
        plt.legend(loc=4, fontsize=18)
        plt.savefig(self.analysispath + '/central_90_dNdE.png',bbox_inches='tight')

        self.energy_range = (np.min([low_5_min_dec, low_5_max_dec]),
                             np.max([high_5_min_dec, high_5_max_dec]))
        print('90% Central Energy Range: {}, {} GeV'.format(round(self.energy_range[0]), round(self.energy_range[1])))
        self.save_items['energy_range'] = self.energy_range

    def ns_scan(self):
        """ NS SCAN WITH SPATIAL PRIOR NOT AN OPTION
        """
        print("ns scan not an option for skymap based analyses")

    def write_circular(self):
        pass

        
class PointSourceFollowup(FastResponseAnalysis):
    """
    Class for point-source or extended source followup
    i.e. there is a fixed location on the sky, not a healpy skymap
    """
    _nside = 256
    def __init__(self, name, ra, dec, tstart, tstop, extension=None,
                 skipped=None, outdir=None, save=True, seed=None):
        
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
        inj.fill(
            self.dec,
            self.llh.exp,
            self.llh.mc,
            self.llh.livetime,
            temporal_model=self.llh.temporal_model)
        self.inj = inj

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
            print("ns = {}\n\n".format(ns))
        self.ts, self.ns = ts, ns
        self.save_items['ts'] = ts
        self.save_items['ns'] = ns
        return ts, ns

    def find_coincident_events(self, ns_params=None):
        r"""Find "coincident events" for the analysis.
        These are ontime events that have a spatial times energy weight greater than 10

        Parameters
        -----------
        ns_params: dict
            Fit parameters to use for weight calculation, e.g. if fit happened in MultiPointSourceFollowup
        """
        if ns_params is None:
            ns_params = self.ns_params
        spatial_weights = self.llh.llh_model.signal(
            self.ra, self.dec, self.llh._events, 
            src_extension=self.extension)[0] / self.llh._events['B']
        energy_ratio, _ = self.llh.llh_model.weight(
            self.llh._events, **ns_params)
        temporal_weights = self.llh.temporal_model.signal(self.llh._events)
        msk = spatial_weights * energy_ratio * temporal_weights > 10
        self.coincident_events = []
        if len(spatial_weights[msk]) > 0:
            self.coincident_events = []
            for ev, s_w, en_w in zip(self.llh._events[msk], 
                        spatial_weights[msk], energy_ratio[msk]):
                self.coincident_events.append({})
                for key in ['run', 'event', 'ra', 'dec', 'sigma', 'logE', 'time']:
                    self.coincident_events[-1][key] = ev[key]
                del_psi = sensitivity_utils.deltaPsi([ev['dec']], [ev['ra']], [self.dec], [self.ra])[0]
                self.coincident_events[-1]['delta_psi'] = del_psi 
                self.coincident_events[-1]['spatial_w'] = s_w
                self.coincident_events[-1]['energy_w'] = en_w
        self.save_items['coincident_events'] = self.coincident_events

    def ns_scan(self, params = {'spectrum': 'dN/dE = 1.00e+00 * (E / 1.00e+03 GeV)^-2.00 [GeV^-1cm^-2s^-1]'}):
        r""" Calculate the llh as a function of ns, to 
        see how shallow the likelihood space is.

        Parameters
        -----------
        params: mappable
            params dictionary for skylab.ps_llh object. Default is 
            that of this analysis
            
        Returns
        ---------
        xs, delta_llh: array-like
            Values of ns scanned and corresponding -2*delta_llh values"""
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
                ax.text(3.5, 0.8, r'Sens. = {:.2f} events'.format(fit_dict['sens']), fontsize = 16)
                ax.text(3.5, 0.7, r' = {:.1e}'.format(self.upperlimit * self.duration * 86400. * 1e6) + r' GeV cm$^{-2}$', fontsize = 16)
                #ax.text(5, 0.5, r'Sens. = {:.2e}'.format(self.inj.mu2flux(fit_dict['sens'])) + ' GeV^-1 cm^-2 s^-1')
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

    def make_dNdE(self):
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
        plt.xlim(1e1, 1e8)
        plt.legend(loc=4, fontsize=18)
        plt.savefig(self.analysispath + '/central_90_dNdE.png',bbox_inches='tight')

        print('90% Central Energy Range: {}, {} GeV'.format(round(low_5), round(high_5)))
        self.save_items['energy_range'] = (self.low5, self.high5)

    def write_circular(self):
        pass

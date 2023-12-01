from multiprocessing import Value

from numpy.lib.recfunctions   import append_fields
import datetime, os
import numpy as np
import healpy as hp
from astropy.time import Time
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import pickle

from .FastResponseAnalysis import PriorFollowup
from .reports import GravitationalWaveReport
from skylab.datasets        import Datasets
from skylab.llh_models      import EnergyLLH
from skylab.spectral_models import PowerLaw 
from skylab.temporal_models import BoxProfile, TemporalModel
from skylab.ps_llh          import PointSourceLLH
import fast_response
from . import sensitivity_utils

class GWFollowup(PriorFollowup):
    '''
    Class for followup of a GW.
    By default, fits the index in the LLH.

    See also:
    ----------
    PriorFollowup: class for skymap-based analyses
    '''

    _dataset = 'GFUOnline_v001p02'
    _fix_index = False
    _float_index = True
    _index = 2.0
    _pixel_scan_nsigma = 3.0
    _allow_neg = True
    _containment = None
    #_nside = 256
    _season_names = ['IC86, 2017', 'IC86, 2018', 'IC86, 2019']
    _nb_days = 5.
    _ncpu = 10

    def __init__(self, name, skymap_path, tstart, tstop, skipped=None, seed=None,
                 outdir=None, save=True, extension=None):
        if (Time(tstop, format='iso').mjd - Time(tstart, format='iso').mjd) < (1001./86400.):
            #Raamis uses 512 for bg trials
            self._nside = 512
        else: 
            #use 256 for 2 week
            self._nside = 256
        super().__init__(name, skymap_path, tstart, tstop, skipped=skipped, seed=seed,
                       outdir=outdir, save=save, extension=extension)
    
    def run_background_trials(self, month=None, ntrials=1000):
        r'''For GW followups with specific time windows,
        just use precomputed background arrays.
        2 allowed time windows:
         - 1000s: for all (default)
         - [-0.1, +14]: BNS, NSBH
        Returns:
        --------
        tsd: array-like
            test-statistic distribution with weighting 
            from alert event spatial prior
        '''

        if self.duration > 1.:
            #self.duration is in days, so this is to load 2 week precomputed trials
            #This is an exact copy of AlertFollowup.py, as it was run with the same script
            from scipy import sparse
            current_rate = self.llh.nbackground / (self.duration * 86400.) * 1000.

            closest_rate = sensitivity_utils.find_nearest(np.linspace(6.0, 7.2, 7), current_rate)
            print(f'Loading 2 week bg trials, rate: {closest_rate}')

            bg_trial_dir = '/data/ana/analyses/NuSources/' \
                + '2023_realtime_gw_analysis/fast_response/' \
                + 'precomputed_background/'

            pre_ts_array = sparse.load_npz(
                bg_trial_dir
                + 'gw_precomputed_trials_delta_t_'
                + '{:.2e}_trials_rate_{:.1f}_low_stats.npz'.format(
                self.duration * 86400., closest_rate, self.duration * 86400.))
            
            #check for nside mismatch
            #if hp.pixelfunc.get_nside(pre_ts_array.shape[1]) != self.nside:
            #    print('Error! Analysis uses nside of %i '.format(self.nside)+
            #          'while precomputed BG is nside %i'.format(hp.pixelfunc.get_nside(pre_ts_array.shape[1])))

            ts_norm = np.log(np.amax(self.skymap))
            ts_prior = pre_ts_array.copy()
            ts_prior.data += 2.*(np.log(self.skymap[pre_ts_array.indices]) - ts_norm)
            ts_prior.data[~np.isfinite(ts_prior.data)] = 0.
            ts_prior.data[ts_prior.data < 0] = 0.
            tsd = ts_prior.max(axis=1).A
            tsd = np.array(tsd)
            self.tsd = tsd
            return tsd

        else: 
            #Case: load 1000s precomputed trials (run by Raamis in O3)
            print('Loading bg trials for 1000s follow-up')
            if month is None:
                # month = datetime.datetime.utcnow().month
                month = Time(self.centertime, format='mjd').datetime.month

            bg_trial_dir = '/data/ana/analyses/NuSources/' \
                + '2021_v2_alert_stacking_FRA/fast_response/gw_precomputed_trials/'

            pre_ts_array = np.load(
                f'{bg_trial_dir}ts_map_{month:02d}.npy',
                allow_pickle=True,
                encoding='latin1'
            )
            if 512 != self.nside:
                print('Error! Analysis uses nside of %i '.format(self.nside)+
                      'while precomputed BG is nside 512')

            # Create spatial prior weighting
            max_ts = []
            ts_norm = np.log(np.amax(self.skymap))
            for i in range(pre_ts_array.size):
                # If a particular scramble in the pre-computed ts_array is empty,
                #that means that sky scan had no events in the sky, so max_ts=0
                if pre_ts_array[i]['ts'].size == 0:
                    max_ts.append(-1*np.inf)
                else:
                    theta, ra = hp.pix2ang(512, pre_ts_array[i]['pixel'])
                    interp = hp.get_interp_val(self.skymap, theta, ra)
                    interp[interp<0] = 0.
                    ts_prior = pre_ts_array[i]['ts'] + 2*(np.log(interp) - ts_norm)
                    max_ts.append(ts_prior.max())
        
            max_ts = np.array(max_ts)

            self.tsd = max_ts
            return max_ts

    def check_events_after(self):
        '''Check that we have recieved at least one event
        after the on-time window closes. 
        This is to make sure we have all data from i3Live before running.

        Returns:
        ----------
        bool: 
         - True if there is at least 1 event after tw (prints how many)
         - False if there are no events after (should re-load!) 
        '''

        t1 = Time(datetime.datetime.utcnow()).mjd
        if ((t1-self.stop)*86400.)>5000.:
            #if it's been long enough, only load 1000s
            print('Loading 2000s of data after the time window')
            t1 = self.stop + 2000./86400. 
        exp_long, livetime_long, grl_long = self.dset.livestream(
            self.start,
            t1,
            load_mc=False,
            floor=self._floor
            )
        
        mask = (exp_long['time']>self.stop)
        
        check_passed = False
        if exp_long[mask].size > 0:
            check_passed = True
            print('Found {} events after end of time window'.format(exp_long[mask].size))
        return check_passed

    def initialize_llh(self, skipped=None, scramble=False):
        '''
        Grab data and format it all into a skylab llh object -
        This function is very similar to the one in FastResponseAnalysis.py, 
        with the main difference being that these are low-latency enough that we may
        have to wait for the last min of data to reach i3live.
        if so, initialize LLH, load ontime data, then add temporal info and ontime data to llh.

        Parameters:
        -----------
        skipped: array of tuples
            Event(s) to be attempted to be removed in the analysis. Format: [(run_id,event_id),... ]
        scramble: bool
            scramble events in LLH (default False)
        
        Returns:
        Skylab LLH object

        See also:
        ----------
        FastResponseAnalysis.initialize_llh: similar function, 
        used throughout the rest of FRA
        get_data: load data function
        check_events_after: looks for events after the ontime window closes
        remove_event: remove a particular event from LLH
        '''
        t0 = Time(datetime.datetime.utcnow()).mjd

        if self.stop > t0 + 60./86400.:
            self.get_data(livestream_start=self.start-6., livestream_stop=self.start)
            print('Loading off-time data')
        elif self.exp is None:
            dset = Datasets[self.dataset]
            self.dset = dset
            check_passed = False
            print('Checking for events after time window')
            while not check_passed:
                check_passed = self.check_events_after()

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
            #timescramble=True,             # not just RA scrambling
            llh_model=llh_model,           # likelihood model
            nsource_bounds=(0., 1e3),      # bounds on fitted ns
            src_extension = src_extension, # Model symmetric extended source
            nsource=1.,                    # seed for nsignal fit
            seed=self.llh_seed)           

        if self.stop > t0 + 60./86400.:
            check_passed = False
            print('Checking for events after time window')
            while not check_passed:
                check_passed = self.check_events_after()

            exp_on, livetime_on, grl_on = self.dset.livestream(
                self.start, 
                self.stop,
                load_mc=False, 
                floor=self._floor,
                wait_until_stop=True)

            llh.append_exp(exp_on,livetime_on)
            self.grl = np.concatenate([self.grl, grl_on])
            self.exp = np.concatenate([self.exp,exp_on])
            self.livetime = self.grl['livetime'].sum()

        box = TemporalModel(
            grl=self.grl,
            poisson_llh=True,
            days=self._nb_days,
            signal=BoxProfile(self.start, self.stop))

        llh.set_temporal_model(box)

        return llh
    
    def get_best_fit_contour(self, proportions=[0.5,0.9]):
        '''Get a contour for the error around the best-fit point. 
        Makes a zoomed skymap of the ts-space with contours

        Parameters:
        -----------
        proportions: list
            confidence levels to calculate (default [0.5,0.9])
        '''
        if self.tsd is None: return
        
        from scipy import special
        from . import plotting_utils
        import meander
        import seaborn as sns
        
        print('Calculating contour around best-fit TS')

        #get threshold TS value for that level in the bg distribution
        dof = 2
        delta_llh_levels = special.gammaincinv(dof/2.0, np.array([(1-prop) for prop in proportions]))
        levels=2*delta_llh_levels

        #for level in levels:
        #    print(level)
        #    print(len(self.ts_scan['TS_spatial_prior_0'][self.ts_scan['TS_spatial_prior_0']>level]))
        with open(self.analysispath + '/fit_ts_map.pkl','wb') as f:
            pickle.dump(self.ts_scan,f)

        #sample_points = np.array(hp.pix2ang(self.nside, np.arange(len(self.skymap)))).T
        loc=np.array((np.pi/2 - self.ts_scan['dec'], self.ts_scan['ra'])).T
        contours_by_level = meander.spherical_contours(loc, self.ts_scan['TS_spatial_prior_0'], levels)

        thetas = []; phis=[]
        for contours in contours_by_level:
            for contour in contours:
                theta, phi = contour.T
                phi[phi<0] += 2.0*np.pi
                thetas.append(theta)
                phis.append(phi)

        #norm_ts = self.ts_scan['TS_spatial_prior_0'] / sum(self.ts_scan['TS_spatial_prior_0'])
        #thetas, phis = plotting_utils.plot_contours(proportions, norm_ts)
        
        #make the plot
        pdf_palette = sns.color_palette("Blues", 500)
        cmap = mpl.colors.ListedColormap(pdf_palette)

        plotting_utils.plot_zoom(self.ts_scan['TS_spatial_prior_0'], self.skymap_fit_ra, self.skymap_fit_dec,
                                 "", range = [0,10], reso=3., cmap = cmap)
        
        plotting_utils.plot_color_bar(labels=[0,round(max(self.ts_scan['TS_spatial_prior_0']))], cmap=cmap, col_label=r"TS",offset=-150)
        cont_ls = ['solid', 'dashed']*(int(len(proportions)/2))
        cont_labels=[f'{proportion*100:.0f}/% CL' for proportion in proportions]

        for i in range(len(thetas)):
            hp.projplot(thetas[i], phis[i], linewidth=2., c='k')#, ls=cont_ls[i], label=cont_labels[i])

        plt.scatter(0,0, marker='*', c = 'k', s = 130, label = "Scan Hot Spot") 
        plt.legend(loc = 2, ncol=1, fontsize = 16, framealpha = 0.95)

        plt.savefig(self.analysispath + '/' + self.analysisid + 'ts_contours.png',bbox_inches='tight')
        plt.savefig(self.analysispath + '/' + self.analysisid + 'ts_contours.pdf',bbox_inches='tight', dpi=300)
        plt.close()
        
    def plot_ontime(self, with_contour=True, contour_files=None, label_events=True):
        return super().plot_ontime(with_contour=True, contour_files=contour_files, label_events=label_events)

    def write_circular(self):
        '''
        Generate a circular from a template. Uses a different template if
        high significance or low-significance (no longer sent for O4).
        High-significance are generated for p<0.2, in case LLAMA sees p<0.01
        for the same event and a circular is needed.
        Saves a text file in the output directory as gcn_[eventname].txt
        '''
        base = os.path.dirname(fast_response.__file__)
        events = self.coincident_events
        pvalue = self.p
        start_iso = str(Time(self.start, format = 'mjd', scale = 'utc').iso)
        stop_iso = str(Time(self.stop, format = 'mjd', scale = 'utc').iso)
        
        namelist = self.name.split('-')
        gw_name = namelist[0]
        try:
            noticeID = namelist[1]+'-'+namelist[2]
        except:
            noticeID = 'NOTICEID'

        if pvalue > 0.2:
            template_path = os.path.join(base, 'circular_templates/gw_gcn_template_low.txt')
        else:
            template_path = os.path.join(base, 'circular_templates/gw_gcn_template_high.txt')

        if pvalue>0.2:
            with open(template_path, 'r') as gcn_template:

                gcn = gcn_template.read()
                low_sens, high_sens = self.sens_range
                gcn = gcn.replace('<lowSens>', '{:1.3f}'.format(low_sens))
                gcn = gcn.replace('<highSens>', '{:1.3f}'.format(high_sens))
                gcn = gcn.replace('<name>' , gw_name)
                gcn = gcn.replace('<tstart>', start_iso)
                gcn = gcn.replace('<tstop>', stop_iso)
                gcn = gcn.replace('<noticeID>', noticeID)

            gcn_file = open(self.analysispath + '/gcn_%s.txt' % gw_name, 'w')
            gcn_file.write(gcn)
            gcn_file.close()

        else:
            #significance = '{:1.2f}'.format(self.significance(pvalue))

            info = ' <dt>\t <ra>\t\t <dec>\t\t <angErr>\t\t\t <p_gwava>\t\t\t <p_llama>\n'
            table = ''
            n_coincident_events=0
            for event in events:
                if event['pvalue']<=0.1:
                    ev_info = info
                    ra = '{:.2f}'.format(np.rad2deg(event['ra']))
                    dec = '{:.2f}'.format(np.rad2deg(event['dec']))
                    sigma = '{:.2f}'.format(np.rad2deg(event['sigma']*2.145966))
                    dt = '{:.2f}'.format((event['time']-self.centertime)*86400.)
                    ev_info = ev_info.replace('<dt>', dt)
                    ev_info = ev_info.replace('<ra>', ra)
                    ev_info = ev_info.replace('<dec>', dec)
                    ev_info = ev_info.replace('<angErr>', sigma)
                    if event['pvalue']<0.0013499:
                        pval_str = '<0.00135'
                        ev_info = ev_info.replace('<p_gwava>', pval_str)
                    else:
                        pval_str = '{:1.3f}'.format(event['pvalue'])
                        ev_info = ev_info.replace('<p_gwava>', pval_str)
                    # ev_info = ev_info.replace('<p_gwava>','{:.3f}'.format(pvalue))
                    table+=ev_info
                    n_coincident_events+=1

            if n_coincident_events == 0:
                num = 'Zero'
            elif n_coincident_events == 1:
                num = 'One'
            elif n_coincident_events == 2:
                num = 'Two'
            else:
                num = str(n_coincident_events)
            #num = events['pvalue'][events['pvalue']<=0.1].size
            gcn_file = open(self.analysispath+'/gcn_%s.txt' % gw_name,'w')
            best_ra = '{:1.2f}'.format(np.rad2deg(self.skymap_fit_ra))
            best_dec = '{:1.2f}'.format(np.rad2deg(self.skymap_fit_dec))

            with open(template_path,'r') as gcn_template:

                for line in gcn_template.readlines():
                    line = line.replace('<N>', num)
                    line = line.replace('<name>', gw_name)
                    line = line.replace('<noticeID>', noticeID)
                    line = line.replace('<tstart>', start_iso)
                    line = line.replace('<tstop>', stop_iso)
                    line = line.replace('<best_ra>', best_ra)
                    line = line.replace('<best_dec>', best_dec)
                    
                    if pvalue<0.0013499:
                        pval_str = '<0.00135'
                        line = line.replace('<p_gwava>', pval_str)
                        #line = line.replace('<sig_gwava>', '>3')
                    else:
                        pval_str = '{:1.3f}'.format(pvalue)
                        line = line.replace('<p_gwava>', pval_str)
                        #line = line.replace('<sig_gwava>', significance)

                    if '<dt>' in line:
                        line = table

                    gcn_file.write(line)
                gcn_file.close()

    def inject_scan(self, ra, dec, ns, poisson=True):
        r''' Run All sky scan using event localization as 
        spatial prior, while also injecting events according
        to event localization

        Parameters:
        -----------
        ns: float
            Number of signal events to inject
        poisson: bool
            Will poisson fluctuate number of signal events
            to be injected

        Returns:
        ----------
        tuple: (results, events)
            results dict with ts, ns, gamma, ra, dec; 
            events in the time window
        '''
        from skylab.priors import SpatialPrior
        from skylab.ps_injector import PointSourceInjector

        ### Set up spatial prior to be used in scan
        spatial_prior = SpatialPrior(self.skymap, allow_neg=self._allow_neg, containment=self._containment)
        pixels = np.arange(len(self.skymap))

        ## Perform all sky scan
        inj = PointSourceInjector(gamma=2, E0=1000.)
        inj.fill(dec, self.llh.exp, self.llh.mc, self.llh.livetime,
                 temporal_model=self.llh.temporal_model)
        ni, sample = inj.sample(ra,ns,poisson=poisson)
        print('injected neutrino at:')
        print(np.rad2deg(sample['ra']),np.rad2deg(sample['dec']))

        val = self.llh.scan(0.0, 0.0, scramble = False, spatial_prior=spatial_prior,
                            time_mask = [self.duration/2., self.centertime],
                            pixel_scan=[self.nside, self._pixel_scan_nsigma],inject=sample)

        exp = self.llh.inject_events(self.llh.exp, sample)
        exp_theta = 0.5*np.pi - exp['dec']
        exp_phi   = exp['ra']
        exp_pix   = hp.ang2pix(self.nside, exp_theta, exp_phi)
        overlap   = np.isin(exp_pix, self.ipix_90)

        t_mask = (exp['time'] <= self.stop) & (exp['time'] >= self.start)
        events = exp[t_mask]

        # add field to see if neutrino is within 90% GW contour
        events = append_fields(events, names=['in_contour', 'ts', 'ns', 'gamma', 'B'],
                              data=np.empty((5, events['ra'].size)),
                              usemask=False)

        for i in range(events['ra'].size):
            events['in_contour'][i]=overlap[i]

        for i in range(events['ra'].size):
            events['B'][i] = self.llh.llh_model.background(events[i])

        if val['TS'].size==0:
            return (0, 0, 2.0, None)
        else:
            ts=val['TS_spatial_prior_0'].max()
            maxLoc = np.argmax(val['TS_spatial_prior_0'])
            ns=val['nsignal'][maxLoc]
            gamma=val['gamma'][maxLoc]
            ra = val['ra'][maxLoc]
            dec = val['dec'][maxLoc]

        val_pix = hp.ang2pix(self.nside,np.pi/2.-val['dec'],val['ra'])
        for i in range(events['ra'].size):
            idx, = np.where(val_pix==exp_pix[t_mask][i])
            events['ts'][i] = val['TS_spatial_prior_0'][idx[0]]
            events['ns'][i] = val['nsignal'][idx[0]]
            events['gamma'][i] = val['gamma'][idx[0]]

        results = dict([('ts',ts),('ns',ns),('gamma',gamma),('ra',ra),('dec',dec)])
        return (results, events)

    def find_coincident_events(self):
        r'''Find "coincident events" for a skymap
        based analysis. These are ALL ontime events,
        with a bool to indicate if they are in the 90% contour

        See Also: 
        ----------
        PriorFollowup.find_coincident_events: 
            similar function for other prior followups
        '''
        if self.ts_scan is None:
            raise ValueError("Need to unblind TS before finding events")
        exp_theta = 0.5*np.pi - self.llh.exp['dec']
        exp_phi   = self.llh.exp['ra']
        exp_pix   = hp.ang2pix(self.nside, exp_theta, exp_phi)

        t_mask=(self.llh.exp['time'] <= self.stop) & (self.llh.exp['time'] >= self.start)
        events = self.llh.exp[t_mask]
        ontime_pix = hp.ang2pix(self.nside, 0.5*np.pi - events['dec'], events['ra'])
        overlap    = np.isin(ontime_pix, self.ipix_90)

        events = append_fields(
            events, names=['in_contour', 'ts', 'ns', 'gamma', 'B'],
            data=np.empty((5, events['ra'].size)),
            usemask=False)

        for i in range(events['ra'].size):
            events['in_contour'][i]=overlap[i]
            events['B'][i] = self.llh.llh_model.background(events[i])

        val_pix = self.scanned_pixels
        for i in range(events['ra'].size):
            idx, = np.where(val_pix == exp_pix[t_mask][i])
            events['ts'][i] = self.ts_scan['TS_spatial_prior_0'][idx[0]]
            events['ns'][i] = self.ts_scan['nsignal'][idx[0]]
            events['gamma'][i] = self.ts_scan['gamma'][idx[0]]

        self.events_rec_array = events
        self.coincident_events = [dict(zip(events.dtype.names, x)) for x  in events]
        self.save_items['coincident_events'] = self.coincident_events

    def upper_limit(self):
        r'''Get a *Sensitivity Range* (not truly an UL)
        for the full range of declinations for GW skymap

        See also:
        ----------
        ps_sens_range: function to get the point source sensitivity
        range within the 90% contour of the map
        make_dec_pdf: plot declination PDF of skymap with IceCube PS sensitivity
        '''
        sens_range = self.ps_sens_range()
        self.sens_range = sens_range
        self.save_items['sens_range'] = sens_range
        self.make_dec_pdf()

    def ps_sens_range(self):
        r''' Compute minimum and maximum sensitivities within
        the declination range of the 90% contour of a given skymap

        Returns:
        --------
        low: float
            lowest sensitivity within dec range
        high: float
            highest sensitivity within dec range
        '''
        
        sens_dir = '/data/ana/analyses/NuSources/2023_realtime_gw_analysis/' \
                +  'fast_response/ps_sensitivities'

        with open(f'{sens_dir}/ps_sensitivities_deltaT_{self.duration*86400.:.2e}s.pkl','rb') as f:
            saved_sens=pickle.load(f)
            dec_range=saved_sens['dec']
            sens=saved_sens['sens_flux']
        #dec_range = np.linspace(-85,85,35)
        #sens = [1.15, 1.06, .997, .917, .867, .802, .745, .662,
        #        .629, .573, .481, .403, .332, .250, .183, .101,
        #        .035, .0286, .0311, .0341, .0361, .0394, .0418,
        #        .0439, .0459, .0499, .0520, .0553, .0567, .0632,
        #        .0679, .0732, .0788, .083, .0866]

        src_theta, src_phi = hp.pix2ang(self.nside, self.ipix_90)
        src_dec = np.pi/2. - src_theta
        src_dec = np.unique(src_dec)

        sens = np.interp(np.degrees(src_dec),dec_range,sens)
        low = sens.min()
        high = sens.max()

        return low,high

    def per_event_pvalue(self):
        '''Calculate per-event p-values. 
        There are a few cases here: 

        * overall p < 0.1: 
        Redoes the all-sky scan, but with only that single event.
        This is the same as asking the quesiton: 
        If that single event is the only one on the sky, with this given skymap,
        what TS/p-value would we get for that event?

        * 1.0 > overall p > 0.1:
        Calculates the p-value at the reconstructed event direction. 
        Takes the TS at that location, and calculates the p-value at that location. 
        Does not re-run the scan, to save time in realtime

        * p=1.0:
        Does not get p-values for the events (all are set to None)

        See also:
        -----------
        per_event_scan: runs the scan, for each individual event
        '''
        self.events_rec_array = append_fields(
            self.events_rec_array,
            names=['pvalue'],
            data=np.empty((1, self.events_rec_array['ra'].size)),
            usemask=False
        )

        if self.p < 0.1:
            for i in range(self.events_rec_array.size):
                ts, p = self.per_event_scan(self.events_rec_array[i])
                self.events_rec_array['pvalue'][i] = p
        elif self.tsd is None:
            for i in range(self.events_rec_array.size):
                self.events_rec_array['pvalue'][i] = None
        else:
            for i in range(self.events_rec_array.size):
                p = np.count_nonzero(self.tsd >= self.events_rec_array['ts'][i]) / float(len(self.tsd))
                self.events_rec_array['pvalue'][i] = p
        
        self.coincident_events = [dict(zip(self.events_rec_array.dtype.names, x)) for x  in self.events_rec_array]
        self.save_items['coincident_events'] = self.coincident_events

    def per_event_scan(self, custom_events):
        '''Runs the all-sky scan for only one (or certain) events on the sky

        Parameters:
        ------------
        custom_events: Skylab event(s), masked array
            Ontime event (or events) to use when running the all sky scan
        
        Returns:
        -----------
        ts, p: float, float
            TS and p-value for the given event(s)
        
        See also:
        ----------
        called by the per_event_pvalue function for individual events on the sky
        '''
        from skylab.priors import SpatialPrior

        spatial_prior = SpatialPrior(self.skymap, containment = self._containment, allow_neg=self._allow_neg)
        val = self.llh.scan(
            0.0,0.0, scramble = False, spatial_prior=spatial_prior,
            time_mask = [self.duration/2., self.centertime],
            pixel_scan=[self.nside, self._pixel_scan_nsigma],
            custom_events=custom_events
        )
        if val['TS'].size == 0:
            ts = -1.*np.inf
        else:
            ts = val['TS_spatial_prior_0'].max()
        p = np.count_nonzero(self.tsd >= ts) / float(len(self.tsd))
        return ts, p

    def make_dec_pdf(self):
        r''' Plot PDF of source declination overlaid with IceCube's
        point source sensitivity
        '''

        sinDec_bins = np.linspace(-1,1,30)
        bin_centers = (sinDec_bins[:-1] + sinDec_bins[1:]) / 2
        
        sens_dir = '/data/ana/analyses/NuSources/2023_realtime_gw_analysis/' \
                +  'fast_response/ps_sensitivities'
        
        with open(f'{sens_dir}/ps_sensitivities_deltaT_{self.duration*86400.:.2e}s.pkl','rb') as f:
            saved_sens=pickle.load(f)
            dec_range=np.sin(saved_sens['dec']*np.pi/180)
            sens=saved_sens['sens_flux']

        #dec_range = np.linspace(-1,1,35)
        #sens = [1.15, 1.06, .997, .917, .867, .802, .745, .662,
        #        .629, .573, .481, .403, .332, .250, .183, .101,
        #        .035, .0286, .0311, .0341, .0361, .0394, .0418,
        #        .0439, .0459, .0499, .0520, .0553, .0567, .0632,
        #        .0679, .0732, .0788, .083, .0866]
        sens = np.array(sens)

        pixels = np.arange(len(self.skymap))
        theta, ra = hp.pix2ang(self.nside, pixels)
        dec = np.pi/2 - theta
        sindec = np.sin(dec)

        pdf = []
        for i in range(sinDec_bins.size-1):
            dec_min = sinDec_bins[i]; dec_max = sinDec_bins[i+1]
            mask = (sindec < dec_max) & (sindec > dec_min)
            pdf.append(self.skymap[pixels[mask]].sum())

        pdf = np.array(pdf)
        pdf /= np.diff(sinDec_bins)

        fig,ax1 = plt.subplots()
        ax1.set_xlabel('sin($\delta$)')
        ax1.set_ylabel('Probability density')
        ax1.set_xlim(-1,1)
        ax1.plot(bin_centers,pdf, color='C0', label='PDF')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        ax2.set_ylabel('E$^2$F (GeVcm$^2$)')  # we already handled the x-label with ax1
        ax2.plot(dec_range, sens, color='C1', label='PS Sensitivity')
        ax2.set_yscale('log')
        ax2.set_xlim(-1,1)
        ax2.tick_params(axis='y')

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.title('%s' % self.analysisid.replace('_', ' '))
        h1, l1 = ax1.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        plt.legend(h1+h2,l1+l2,loc=1)
        if self.save_output:
            plt.savefig(
                self.analysispath + f'/decPDF.png',
                bbox_inches='tight', dpi=200
            )
        plt.close()

    def generate_report(self):
        r'''Generates GW report using class attributes
        and the ReportGenerator Class
        
        See also:
        ----------
        fast_response/reports: GravitationalWaveReport
        '''
        report = GravitationalWaveReport(self)
        #if self.duration > 1.:
        #    report._figure_list = [('decPDF', 'decPDF.png'),('ts_contours', 'ts_contours.png')]

        report.generate_report()
        report.make_pdf()
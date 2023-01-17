from multiprocessing import Value

from numpy.lib.recfunctions   import append_fields
import datetime, os
import numpy as np
import healpy as hp
from astropy.time import Time
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from .FastResponseAnalysis import PriorFollowup
from .reports import GravitationalWaveReport
from skylab.llh_models      import EnergyLLH
from skylab.spectral_models import PowerLaw 
from skylab.temporal_models import BoxProfile, TemporalModel
from skylab.ps_llh          import PointSourceLLH
import fast_response
from . import sensitivity_utils

class GWFollowup(PriorFollowup):
    _dataset = 'GFUOnline_v001p02'
    _fix_index = False
    _float_index = True
    _index = 2.0
    _pixel_scan_nsigma = 3.0
    _allow_neg = True
    _containment = None
    _nside = 256
    _season_names = ['IC86, 2017', 'IC86, 2018', 'IC86, 2019']
    _nb_days = 5.
    _ncpu = 10
    
    def run_background_trials(self, month=None, ntrials=1000):
        '''For GW followup, 2 possible time windows:
        1000s: for all (default)
        [-0.1, +14]: NS included (run manually)
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

            #TODO: replace here once done
            #closest_rate = sensitivity_utils.find_nearest(np.linspace(6.0, 7.2, 7), current_rate)
            rates = [6.0,6.2,6.4,6.8,7.0]
            closest_rate = sensitivity_utils.find_nearest(rates, current_rate)
            print(f'Loading 2 week bg trials, rate: {closest_rate}')

            #bg_trial_dir = '/data/ana/analyses/NuSources/' \
            #    + '2023_realtime_gw_analysis/fast_response/' \
            #    + 'precomputed_background/'
            #TODO: change to permanent storage once assigned
            bg_trial_dir = '/data/user/jthwaites/FastResponseAnalysis/output/trials/glob_trials/'
            pre_ts_array = sparse.load_npz(
                bg_trial_dir
                + 'gw_precomputed_trials_delta_t_'
                + '{:.2e}_trials_rate_{:.1f}_low_stats.npz'.format(
                self.duration * 86400., closest_rate, self.duration * 86400.))
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

    def initialize_llh(self, skipped=None, scramble=False):
        '''
        Grab data and format it all into a skylab llh object
        This function is very similar to the one in FastResponseAnalysis.py, 
        with the main difference being that these are low-latency enough that we may
        have to wait for the last min of data to reach i3live.
        Initialize LLH, load ontime data, then add temporal info and ontime data to llh
        '''
        t0 = Time(datetime.datetime.utcnow()).mjd

        if self.stop > t0 + 60./86400.:
            self.get_data(livestream_start=self.start-6., livestream_stop=self.start)
            print('Loading off-time data')
        elif self.exp is None:
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
            exp_on, livetime_on, grl_on = self.dset.livestream(
                self.start, 
                self.stop+60./86400.,
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

    def plot_ontime(self, with_contour=True, contour_files=None):
        return super().plot_ontime(with_contour=True, contour_files=contour_files)

    def write_circular(self):
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

        if pvalue > 0.01:
            template_path = os.path.join(base, 'circular_templates/gw_gcn_template_low.txt')
        else:
            template_path = os.path.join(base, 'circular_templates/gw_gcn_template_high.txt')

        if pvalue>0.01:
            with open(template_path, 'r') as gcn_template:

                gcn = gcn_template.read()
                low_sens, high_sens = self.sens_range
                # for key, val in [('<lowSens>', f'{low_sens:1.3f}'),
                #                  ()]
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
            significance = '{:1.2f}'.format(self.significance(pvalue))

            info = ' <dt>   <ra>       <dec>          <angErr>                    <p_gwava>                 <p_llama>\n'
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
            with open(template_path,'r') as gcn_template:

                for line in gcn_template.readlines():
                    line = line.replace('<N>', num)
                    line = line.replace('<name>', gw_name)
                    line = line.replace('<noticeID>', noticeID)
                    line = line.replace('<tstart>', start_iso)
                    line = line.replace('<tstop>', stop_iso)
                    if pvalue<0.0013499:
                        pval_str = '<0.00135'
                        line = line.replace('<p_gwava>', pval_str)
                        line = line.replace('<sig_gwava>', '>3')
                    else:
                        pval_str = '{:1.3f}'.format(pvalue)
                        line = line.replace('<p_gwava>', pval_str)
                        line = line.replace('<sig_gwava>', significance)

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
        --------
        ts: array
            array of ts values of 
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
        if self.ts_scan is None:
            raise ValueError("Need to unblind TS before finding events")
        exp_theta = 0.5*np.pi - self.llh.exp['dec']
        exp_phi   = self.llh.exp['ra']
        exp_pix   = hp.ang2pix(self.nside, exp_theta, exp_phi)
        overlap   = np.isin(exp_pix, self.ipix_90)

        t_mask=(self.llh.exp['time'] <= self.stop) & (self.llh.exp['time'] >= self.start)
        events = self.llh.exp[t_mask]

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
        high: floaot
            highest sensitivity wihtin dec range
        '''
        dec_range = np.linspace(-85,85,35)
        sens = [1.15, 1.06, .997, .917, .867, .802, .745, .662,
                .629, .573, .481, .403, .332, .250, .183, .101,
                .035, .0286, .0311, .0341, .0361, .0394, .0418,
                .0439, .0459, .0499, .0520, .0553, .0567, .0632,
                .0679, .0732, .0788, .083, .0866]

        src_theta, src_phi = hp.pix2ang(self.nside, self.ipix_90)
        src_dec = np.pi/2. - src_theta
        src_dec = np.unique(src_dec)

        sens = np.interp(np.degrees(src_dec),dec_range,sens)
        low = sens.min()
        high = sens.max()

        return low,high

    def per_event_pvalue(self):
        self.events_rec_array = append_fields(
            self.events_rec_array,
            names=['pvalue'],
            data=np.empty((1, self.events_rec_array['ra'].size)),
            usemask=False
        )

        if self.p < 0.05:
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

        dec_range = np.linspace(-1,1,35)
        sens = [1.15, 1.06, .997, .917, .867, .802, .745, .662,
                .629, .573, .481, .403, .332, .250, .183, .101,
                .035, .0286, .0311, .0341, .0361, .0394, .0418,
                .0439, .0459, .0499, .0520, .0553, .0567, .0632,
                .0679, .0732, .0788, .083, .0866]
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
        r'''Generates report using class attributes
        and the ReportGenerator Class'''
        report = GravitationalWaveReport(self)
        report.generate_report()
        report.make_pdf()
from multiprocessing import Value
from .FastResponseAnalysis import PriorFollowup

from numpy.lib.recfunctions   import append_fields
import datetime
import numpy as np
import healpy as hp
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from .reports import GravitationalWaveReport

class GWFollowup(PriorFollowup):
    _dataset = 'GFUOnline_v001p02'
    _fix_index = False
    _float_index = True
    _index = 2.0
    _pixel_scan_nsigma = 3.0
    _allow_neg = True
    _containment = None
    _nside = 512
    _season_names = ['IC86, 2017', 'IC86, 2018', 'IC86, 2019']
    _nb_days = 5.
    
    def run_background_trials(self, month=None, ntrials=1000):
        if month is None:
            month = datetime.datetime.utcnow().month

        pre_ts_array = np.load(
            f'/data/user/rhussain/ligo_skymaps/ts_map_{month:02d}.npy',
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

    def plot_ontime(self, with_contour=True, contour_files=None):
        return super().plot_ontime(with_contour=True, contour_files=contour_files)

    def write_circular(self):
        pass

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
        else:
            for i in range(self.events_rec_array.size):
                p = np.count_nonzero(self.tsd >= self.events_rec_array['ts'][i]) / float(len(self.tsd))
                self.events_rec_array['pvalue'][i] = p
        
        self.coincident_events = [dict(zip(self.events_rec_array.dtype.names, x)) for x  in self.events_rec_array]
        self.save_items['coincident_events'] = self.coincident_events

    def per_event_scan(self, custom_events):
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
        theta, ra = hp.pix2ang(self.nside,pixels)
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
        ax1.set_ylabel('Probability')
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
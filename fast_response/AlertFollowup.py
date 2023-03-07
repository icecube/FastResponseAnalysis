import numpy as np
import healpy as hp
from scipy import sparse
import pickle, os
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.time import Time

import fast_response
from .FastResponseAnalysis import PriorFollowup
from .reports import AlertReport
from . import sensitivity_utils

class AlertFollowup(PriorFollowup):
    _smear = False
    _dataset = "GFUOnline_v001p02"
    _fix_index = True
    _float_index = not _fix_index
    _index = 2.5

    def run_background_trials(self, ntrials = 1000):
        r'''For alert events with specific time windows,
        just use precomputed background arrays
        
        Returns:
        --------
            tsd: array-like
                test-statistic distribution with weighting 
                from alert event spatial prior
        '''
        current_rate = self.llh.nbackground / (self.duration * 86400.) * 1000.
        closest_rate = sensitivity_utils.find_nearest(np.linspace(6.2, 7.2, 6), current_rate)

        #if self.on_cobalt:
        bg_trial_dir = '/data/ana/analyses/NuSources/' \
            + '2021_v2_alert_stacking_FRA/fast_response/' \
            + 'alert_precomputed_trials/'
        #else:
        #bg_trial_dir = '/cvmfs/icecube.opensciencegrid.org/users/jthwaites/' \
        #   + '2021_v2_alert_stacking_FRA/fast_response/' \
        #    + 'alert_precomputed_trials/'
        
        pre_ts_array = sparse.load_npz(
            bg_trial_dir
            + 'precomputed_trials_delta_t_'
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

    def upper_limit(self):
        sens_range = self.ps_sens_range()
        self.sens_range = sens_range
        self.save_items['sens_range'] = sens_range
        self.sens_range_plot()

    def ps_sens_range(self):
        r'''Compute the minimum and maximum sensitivities
        within the 90% contour of the skymap'''

        #if self.on_cobalt:
        sens_dir = '/data/ana/analyses/NuSources/2021_v2_alert_stacking_FRA/' \
            + 'fast_response/reference_sensitivity_curves/'
        #else:
        #sens_dir = '/cvmfs/icecube.opensciencegrid.org/users/jthwaites/'\
        #   + 'fast_response/reference_sensitivity_curves/'

        with open(f'{sens_dir}ideal_ps_sensitivity_deltaT_{self.duration:.2e}_50CL.pkl', 'rb') as f:
            ideal = pickle.load(f, encoding='bytes')
        delta_t = self.duration * 86400.
        src_theta, src_phi = hp.pix2ang(self.nside, self.ipix_90)
        src_dec = np.pi/2. - src_theta
        src_dec = np.unique(src_dec)
        src_dec = np.sin(src_dec)
        sens = np.interp(src_dec, ideal[b'sinDec'], np.array(ideal[b'sensitivity'])*delta_t*1e6)
        low = sens.min(); high=sens.max()
        return low, high

    def sens_range_plot(self):
        r'''For alert events, make a sensitivity plot highlighting
        the region where the contour lives'''
        fig, ax = plt.subplots()

        #if self.on_cobalt:
        sens_dir = '/data/ana/analyses/NuSources/2021_v2_alert_stacking_FRA/' \
            + 'fast_response/reference_sensitivity_curves/'
        #else:
        #sens_dir = '/cvmfs/icecube.opensciencegrid.org/users/jthwaites/' \
        #   +'fast_response/reference_sensitivity_curves/'

        with open(f'{sens_dir}ideal_ps_sensitivity_deltaT_{self.duration:.2e}_50CL.pkl', 'rb') as f:
            ideal = pickle.load(f, encoding='bytes')
        delta_t = self.duration * 86400.
        plt.plot(ideal[b'sinDec'], np.array(ideal[b'sensitivity'])*delta_t*1e6, lw=3, ls='-', 
                 color=sns.xkcd_rgb['dark navy blue'], label = 'P.S. sensitivity', alpha=0.7)
        #FILL BETWEEN RANGE
        src_theta, src_phi = hp.pix2ang(self.nside, self.ipix_90)
        src_dec = np.pi/2. - src_theta
        src_dec = np.unique(src_dec)
        src_dec = np.sin(src_dec)
        ax.axvspan(src_dec.min(), src_dec.max(), alpha=0.3, color=sns.xkcd_rgb['light navy blue'],
                label='90\% contour region')
        plt.text(0.05, 3e1, 'Min sens.: {:.1e}'.format(self.sens_range[0]) + r' GeV cm$^{-2}$')
        plt.text(0.05, 1.5e1, 'Max sens.: {:.1e}'.format(self.sens_range[1]) + r' GeV cm$^{-2}$')
        plt.grid(which='both', alpha=0.2, zorder=1)
        plt.yscale('log')
        plt.xlabel(r'$\sin \delta$', fontsize=20)
        plt.legend(loc=3, frameon=False)
        plt.ylabel(r'$E^2 \frac{dN_{\nu+\bar{\nu}}}{dEdA} \Big|_{\mathrm{1 TeV}}$ (GeV cm$^{-2}$)', fontsize=20)
        plt.ylim(3e-2, 3e2)
        if self.save_output:
            plt.savefig(self.analysispath + '/upper_limit_distribution.png', bbox_inches='tight', dpi=200)

    def generate_report(self):
        r'''Generates report using class attributes
        and the ReportGenerator Class'''
        report = AlertReport(self)
        report.generate_report()
        report.make_pdf()

    def write_circular(self, alert_results):
        base = os.path.dirname(fast_response.__file__)
        template_path = os.path.join(base, 'circular_templates/internal_followup.txt')

        new_f = []
        high_sig = False
        for window in alert_results.keys():
            if alert_results[window]['p'] < 0.01:
                high_sig = True
        
        analysis_1000 = alert_results[1000.]
        analysis_2day = alert_results[172800.]
        
        if 'Cascade' in analysis_1000['name']:
            alert_id=analysis_1000['name'][:23]
        else:
            alert_id = analysis_1000['name'][:15]
        analysis_1000['gcn_num']=str(analysis_1000['gcn_num'])

        if len(analysis_1000['gcn_num'])>6:
            gcn_link='notices_amon_icecube_cascade/'+analysis_1000['gcn_num']+'.amon'
        else:
            gcn_link= 'gcn3/'+ analysis_1000['gcn_num'] +'.gcn3'
        
        if 'coincident_events' not in analysis_1000.keys():
            analysis_1000['coincident_events'] = []
        ev_is_are = 'event is' if len(analysis_1000['coincident_events']) == 1 else 'events are'
        if not high_sig:
            if len(analysis_1000['coincident_events']) == 0:
                coinc_and_p = ''
            elif len(analysis_1000['coincident_events']) == 1:
                coinc_and_p = 'We find that this additional event is well described by atmospheric\n' \
                    + 'background expectations, with a p-value of {:.2f}. '.format(analysis_1000['p'])
            else:
                coinc_and_p = 'We find that these additional {} events are well described by atmospheric\n' \
                    + 'background expectations, with a p-value of {:.2f}. '.format(len(analysis_1000['coincident_events']), analysis_1000['p'])
            long_p_and_lim = 'In this case, we report a p-value of {:.2f},'.format(analysis_2day['p']) \
                    + ' consistent with no significant \nexcess of track events.' 
        else:
            coinc_and_p = 'We accordingly derive a p-value of {:.3f}.'.format(analysis_1000['p'])
            if analysis_1000['p'] < 0.01:
                coinc_and_p = coinc_and_p + ' Due to the coincidences identified in this search, ' \
                    + 'we strongly encourage followup observations.'
            else:
                pass
            if analysis_2day['p'] < 0.01:
                long_p_and_lim = 'In this case, we report a p-value of {:.3f}.'.format(analysis_2day['p']) \
                    + ' Due to the coincidences identified in this search, we strongly encourage followup observations.'
            else:
                long_p_and_lim = 'In this case, we report a p-value of {:.2f},'.format(analysis_2day['p']) \
                    + ' consistent with no significant \nexcess of track events. '
        
        if len(analysis_1000['coincident_events']) == 0:
            nevents = 'zero'
        elif len(analysis_1000['coincident_events']) == 1:
            nevents = 'one'
        elif len(analysis_1000['coincident_events']) == 2:
            nevents = 'two'
        else:
            nevents = str(len(analysis_1000['coincident_events']))

        if '{:.1e}'.format(analysis_1000['sens_range'][0]) == '{:.1e}'.format(analysis_1000['sens_range'][1]):
            short_sens_range = f'is {analysis_1000["sens_range"][0]:.1e}'
        else:
            short_sens_range = f'ranges from {analysis_1000["sens_range"][0]:.1e} to {analysis_1000["sens_range"][1]:.1e}'
        if '{:.1e}'.format(analysis_2day['sens_range'][0]) == '{:.1e}'.format(analysis_2day['sens_range'][1]):
            long_sens_range = f'is {analysis_2day["sens_range"][0]:.1e}'
        else:
            long_sens_range = f'ranges from {analysis_2day["sens_range"][0]:.1e} to {analysis_2day["sens_range"][1]:.1e}'

        keypairs = [('alert_name', alert_id), ('gcn_number', gcn_link), 
                    ('start_utc', Time(analysis_1000['start'], format='mjd').iso), 
                    ('stop_utc', Time(analysis_1000['stop'], format='mjd').iso), 
                    ('long_start_utc', Time(analysis_2day['start'], format='mjd').iso), 
                    ('long_stop_utc', Time(analysis_2day['stop'], format='mjd').iso),
                    ('n_events', nevents), 
                    ('events_is_are', ev_is_are),
                    ('coincident_events_and_p_str', coinc_and_p),
                    ('long_p_and_lim', long_p_and_lim),
                    ('short_sens_range', short_sens_range),
                    ('long_sens_range', long_sens_range),
                    ('low_en', analysis_1000['energy_range'][0]), 
                    ('high_en', analysis_1000['energy_range'][1])]

        with open(template_path, 'r') as f:
            for line in f.readlines():
                for k, r in keypairs:
                    if k in line:
                        if type(r) == str:
                            form_r = r
                        elif type(r) == int:
                            form_r = '{}'.format(r)
                        elif k in ['low_en', 'high_en']:
                            form_r = '{:.0e}'.format(r)
                        elif 'sens' in k:
                            form_r = '{:.1e}'.format(r)
                        else:
                            form_r = '{}'.format(r)
                        line = line.replace('<'+k+'>', form_r)
                new_f.append(line)
        
        with open(self.analysispath + f'/{alert_id}_circular.txt', 'w') as f:
            for line in new_f:
                f.write(line)
        

class TrackFollowup(AlertFollowup):
    _smear = True
    _dataset = "GFUOnline_v001p02"
    _fix_index = True
    _float_index = not _fix_index
    _index = 2.5

    def format_skymap(self, skymap):
        skymap = self.convert_llh_to_prob(skymap)
        skymap = super().format_skymap(skymap)
        return skymap

    def ipixs_in_percentage(self, percentage):
        """Finding ipix indices confined in a given percentage.
        """
        skymap = self.skymap_llh
        if hp.pixelfunc.get_nside(skymap) != self._nside:
            skymap = hp.pixelfunc.ud_grade(skymap, self._nside, pess=True)
        indices = np.r_[:len(skymap)]
        if percentage == 0.9:
            msk = (skymap < 64.2) * (skymap > 0.)
        elif percentage == 0.5:
            msk = (skymap < 22.2) * (skymap > 0.)
        else:
            raise ValueError('Must use 50\% or 90\% containment for alert events')
        msk *= ~np.isnan(skymap)
        msk *= ~np.isinf(skymap)
        ipix = np.asarray(indices[msk], dtype=int)
        return ipix

    def convert_llh_to_prob(self, skymap_fits):
        '''
        This takes a millipede map and converts from the likelihood
        values to a pdf assuming the order-preserving mapping
        we use to account for systematics
        '''
        skymap_llh = skymap_fits.copy()
        self.skymap_llh = skymap_llh
        skymap_fits = np.exp(-1. * skymap_fits / 2.) #Convert from 2LLH to unnormalized probability
        skymap_fits = np.where(skymap_fits > 1e-12, skymap_fits, 0.0)
        skymap_fits = skymap_fits / np.sum(skymap_fits)
        skymap_fits = skymap_fits/skymap_fits.sum()
        if self._smear:
            ninety_msk = skymap_llh < 64.2
            nside = hp.get_nside(skymap_llh)
            cdf = np.cumsum(np.sort(skymap_fits[ninety_msk][::-1]))
            pixs_above_ninety = np.count_nonzero(cdf> 0.1)
            original_ninety_area = hp.nside2pixarea(nside) * pixs_above_ninety
            new_ninety_area = hp.nside2pixarea(nside) * np.count_nonzero(skymap_fits[ninety_msk])
            original_ninety_radius = np.sqrt(original_ninety_area / np.pi)
            new_ninety_radius = np.sqrt(new_ninety_area / np.pi)
            scaled_probs = self._scale_2d_gauss(skymap_fits, original_ninety_radius, new_ninety_radius)
            skymap_fits = scaled_probs
        return skymap_fits

    def _scale_2d_gauss(self, arr, sigma_arr, new_sigma):
        '''
        Helper function for scaling the likelihood space
        according to systematic resimulations
        '''
        tmp = arr**(sigma_arr**2. / new_sigma**2.)/(np.sqrt(2.*np.pi)*new_sigma)* \
                        np.power(np.sqrt(2.*np.pi)*sigma_arr, (sigma_arr**2. / new_sigma**2.)) 
        return tmp / np.sum(tmp)

class CascadeFollowup(AlertFollowup):
    _smear = False
    _dataset = "GFUOnline_v001p02"
    _fix_index = True
    _float_index = not _fix_index
    _index = 2.5
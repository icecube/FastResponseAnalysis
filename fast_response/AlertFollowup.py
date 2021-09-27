from .FastResponseAnalysis import PriorFollowup

class AlertFollowup(PriorFollowup):
    _smear = False

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
        closest_rate = find_nearest(np.linspace(6.2, 7.2, 6), current_rate)
        pre_ts_array = sparse.load_npz(
            '/data/user/apizzuto/fast_response_skylab/fast-response/'
            + 'fast_response/precomputed_background/glob_trials/'
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
        self.sens_range_plot()

    def ps_sens_range(self):
        r'''Compute the minimum and maximum sensitivities
        within the 90% contour of the skymap'''
        if self.alert_event:
            with open('/data/user/apizzuto/fast_response_skylab/dump/ideal_ps_sensitivity_deltaT_{:.2e}_50CL.pkl'.format(self.duration), 'rb') as f:
                ideal = pickle.load(f, encoding='bytes')
            delta_t = self.duration * 86400.
            src_theta, src_phi = hp.pix2ang(self.nside, self.ipix_90)
            src_dec = np.pi/2. - src_theta
            src_dec = np.unique(src_dec)
            src_dec = np.sin(src_dec)
            sens = np.interp(src_dec, ideal[b'sinDec'], np.array(ideal[b'sensitivity'])*delta_t*1e6)
            low = sens.min(); high=sens.max()
            return low, high
        else:
            return None, None

    def sens_range_plot(self):
        r'''For alert events, make a sensitivity plot highlighting
        the region where the contour lives'''
        fig, ax = plt.subplots()
        with open('/data/user/apizzuto/fast_response_skylab/dump/ideal_ps_sensitivity_deltaT_{:.2e}_50CL.pkl'.format(self.duration), 'rb') as f:
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
            plt.s

class TrackFollowup(AlertFollowup):
    _smear = True

    def format_skymap(self, skymap):
        self.convert_llh_to_prob(skymap)
        skymap = super().format_skymap(skymap)
        return skymap

    def convert_llh_to_prob(self, skymap_fits):
        skymap_llh = skymap_fits.copy()
        skymap_fits = np.exp(-1. * skymap_fits / 2.) #Convert from 2LLH to unnormalized probability
        skymap_fits = np.where(skymap_fits > 1e-20, skymap_fits, 0.0)
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
            scaled_probs = scale_2d_gauss(skymap_fits, original_ninety_radius, new_ninety_radius)
            skymap_fits = scaled_probs
        self.skymap = skymap_fits

    def _scale_2d_gauss(self, arr, sigma_arr, new_sigma):
        tmp = arr**(sigma_arr**2. / new_sigma**2.)/(np.sqrt(2.*np.pi)*new_sigma)* \
                        np.power(np.sqrt(2.*np.pi)*sigma_arr, (sigma_arr**2. / new_sigma**2.)) 
        return tmp / np.sum(tmp)

class CascadeFollowup(AlertFollowup):
    _smear = False
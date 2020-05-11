import numpy as np
from glob import glob
import pandas as pd
import pickle, csv, ast, sys
sys.path.append('/data/user/apizzuto/fast_response_skylab/alert_event_followup/FIRESONG/')
from Firesong import firesong_simulation

data_path = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/FIRESONG/Results/'
eff_area_path = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/effective_areas_alerts/'

bg_rates = {'HESE_gold': 0.4, 'HESE_bronze': 0.9, 'GFU_gold': 5.7, 'GFU_bronze': 13.8}

class Universe():
    r'''Given a set of cosmological parameters, calculate neutrino sources in 
    the universe, and find which of those will initiate alerts in IceCube'''
    def __init__():
        pass

    def find_alerts(self):
        r'''Compile background and signal alerts'''
        if self.uni_header is None:
            self.create_universe()
        signal_alerts = self.find_signal_alerts()
        background_alerts = self.find_background_alerts()
        alerts = {'signal': signal_alerts,
                    'background': background_alerts}
        self.alerts = alerts

    def universe_firesong(self):
        pass

    def create_universe(self):
        r'''
        Simulate universe using FIRESONG
        '''
        pass

    def find_background_alerts(self):
        r'''Sample from expected number of background alerts'''
        bg_alerts = {}
        for sample in ['HESE', 'GFU']:
                for cut, level in [('gold', 'tight'), ('bronze', 'loose')]:
                    N = self.rng.poisson(lam = bg_rates[sample + '_' + cut]*self.data_years)
                    sigs = sample_signalness(cut=level, stream='background', size=N)
                    decs = sample_declination(cut=cut, size=N)
                    skymap_inds, skymap_decs = self.sample_skymap(decs)
                    bg_alerts[sample + '_' + cut] = (N, sigs, decs, skymap_inds, skymap_decs)
        self.bg_alerts = bg_alerts
        #Find background alert skymaps
        return bg_alerts

    def dec_band_from_decs(self):
        r'''Alert stream effective areas given only in three declination
        bands, conversion from dec to band done here'''
        bands = np.digitize(self.sources['dec'], bins = [-90., -5., 30., 90.]) - 1.
        self.sources['dec_bands'] = bands

    def N_per_dec_band(self):
        r'''Assuming a flux normalization of 1, calculate mean number of 
        expected events assuming a given spectrum'''
        ens = np.logspace(1., 8., 501)
        tmp = {}
        for sample in ['HESE', 'GFU']:
            for level in ['gold', 'bronze']:
                tmp[sample + '_' + level] = []
                for dec in [-45., 0., 45]:
                    #test_fl = 1.0 #Is it just the numu contribution for HESE?
                    test_fl = 3.0 if sample == 'HESE' else 1.0 #eff A for HESE is all flavor
                    tmp[sample + '_' + level].append(calc_mean_n(ens, dec, stream=sample, 
                        level=level, gamma=self.diffuse_flux_ind*-1., flux_norm=test_fl))
        self.n_per_dec = tmp

    def find_signal_alerts(self):
        r'''With distribution of sources, calculate where the alert events
        are coming from'''
        sig_alerts = {}
        for stream in ['GFU', 'HESE']:
            for cut, lev in [('gold', 'tight'), ('bronze', 'loose')]:
                sig_alerts[stream + '_' + cut] = [(0,0.,0.)]*len(self.sources['dec'])
        #sig_alerts = [(0,0.,'','')]*len(self.sources['dec'])
        for jjj, (src_dec, src_flux, src_bnd) in enumerate(zip(self.sources['dec'], self.sources['flux'], self.sources['dec_bands'])):
            for stream in ['GFU', 'HESE']:
                for cut, lev in [('gold', 'tight'), ('bronze', 'loose')]:
                    nexp = self.n_per_dec[stream + '_' + cut][int(src_bnd)] * src_flux
                    N = self.rng.poisson(lam=nexp)
                    #if N != 0.0:
                    sigs = sample_signalness(cut=lev, stream='signal', size = N) if N != 0 else 0.
                    sig_alerts[stream + '_' + cut][jjj] = (N, sigs, nexp)
        self.sig_alerts = sig_alerts
        return sig_alerts

    def sample_skymap(self, decs):
        r'''Only use real alert event skymap locations'''
        ###########################################################################
        ################ FIX WHEN STEADY TRIALS FINISH #######################################
        ################ shouldn't stop at 19 index in next line
        #####################################################
        map_decs = np.load('/data/user/apizzuto/fast_response_skylab/alert_event_followup/effective_areas_alerts/decs_by_ind.npy')[1][:19]
        sample_decs, idxs = [], []
        if isinstance(decs, float):
            decs = [decs]
        for dec in decs:
            diffs = np.abs(np.sin(map_decs)-np.sin(dec))
            if np.min(diffs) > 0.1:
                idx = find_nearest_ind(map_decs, dec)
                sample_dec = map_decs[idx]
            else:
                nearby_inds = np.argwhere(diffs < 0.1).flatten()
                idx = self.rng.choice(nearby_inds)
                sample_dec = map_decs[idx]
            sample_decs.append(sample_dec)
            idxs.append(idx)
        return sample_decs, idxs

    def find_alert_skymaps(self):
        r'''Iterate over alert locations, find corresponding alert'''
        skymaps = {}
        for stream in self.sig_alerts.keys():
            skymaps[stream] = [None] * len(self.sources['dec'])
            for jjj, (src_dec, src_flux) in enumerate(zip(self.sources['dec'], self.sources['flux'])):
                if self.sig_alerts[stream][jjj][0] == 0:
                    continue
                else:
                    tmp = self.sample_skymap(np.radians(src_dec))
                    skymaps[stream][jjj] = (tmp[0][0], tmp[1][0])
        self.skymaps = skymaps
        return 
          
    def additional_signal_events(self):
        r'''After finding signal events, calculate any events that could
        accompany the alert in a higher effective area, less pure sample'''
        en_bins = np.logspace(2., 9., 501)
        ens = centers(en_bins)
        extra_events = {}
        for stream in ['GFU', 'HESE']:
            for cut, lev in [('gold', 'tight'), ('bronze', 'loose')]:
                extra_events[stream + '_' + cut] = np.zeros(len(self.sources['dec']))
        if self.sig_alerts is None:
            self.find_signal_alerts()
        with open('/data/user/apizzuto/fast_response_skylab/alert_event_followup/' + \
                'effective_areas_alerts/gfu_online_effective_area_spline.npy', 'r') as f:
            gfu_eff_spline = pickle.load(f)
        for stream in self.sig_alerts.keys():
            for jjj, (src_dec, src_flux) in enumerate(zip(self.sources['dec'], self.sources['flux'])):
                if self.sig_alerts[stream][jjj][0] == 0:
                    continue
                else:
                    mu_extra = np.sum(gfu_eff_spline(np.log10(ens), np.sin(src_dec*np.pi / 180)) * \
                     spectrum(ens, gamma = -1.*self.diffuse_flux_ind, flux_norm=src_flux)*np.diff(en_bins)*1e4)
                    N_extra = self.rng.poisson(mu_extra)
                    extra_events[stream][jjj] = N_extra
        self.extra_events = extra_events
        return extra_events

    

class SteadyUniverse(Universe):
    r'''Universe inherited class for steady neutrino sources'''

    def __init__(self, lumi, evol, density, diffuse_flux_norm, diffuse_flux_ind,
                    **kwargs):

        self.evolution = evol
        self.density = density
        self.lumi = lumi
        self.diffuse_flux_norm = diffuse_flux_norm
        self.diffuse_flux_ind = diffuse_flux_ind
        self.alerts, self.sources, self.uni_header = None, None, None
        self.sim_flux = None
        self.data_years = kwargs.pop('data_years', 1.)
        if self.diffuse_flux_norm is not None:
            self.N_per_dec_band() #initializes effective area
        self.manual_lumi = kwargs.pop('manual_lumi', 0.0)
        self.seed = kwargs.pop('seed', None)
        self.rng = np.random.RandomState(self.seed)
        self.timescale = self.data_years * 365. * 86400.

    def universe_firesong(self):
        return firesong_simulation('', density=self.density, Evolution=self.evolution,
                    Transient=False, fluxnorm = self.diffuse_flux_norm,
                    index=self.diffuse_flux_ind, LF = self.lumi, luminosity=self.manual_lumi, seed=self.seed)

    def create_universe(self):
        r'''
        Run FIRESONG for every year of data
        '''
        tmp_fls, tmp_dec, tmp_zs, tmp_tot = [],[],[],0.
        uni = self.universe_firesong()
        tmp_dec.extend(uni['sources']['dec']), tmp_fls.extend(uni['sources']['flux'])
        tmp_zs.extend(uni['sources']['z'])
        tmp_tot += uni['total_flux']
        #fluxes are E^2 dN/dE at 100 TeV, convert now to dN/dE * DeltaT at 1 GeV
        tmp_fls = np.array(tmp_fls)
        #Don't need the 1+z for time-integrated bursts
        tmp_fls *= self.timescale * np.power(1e5, self.diffuse_flux_ind - 2.)
        self.sources = {'dec': np.array(tmp_dec), 'flux': tmp_fls, 'z': np.array(tmp_zs)} 
        self.uni_header = uni['header']
        self.sim_flux = tmp_tot
        self.dec_band_from_decs()


class TransientUniverse(Universe):
    r'''Universe inherited class for short timescale sources'''

    def __init__(self, lumi, evol, density, diffuse_flux_norm, diffuse_flux_ind,
                    **kwargs):

        self.evolution = evol
        self.density = density
        self.lumi = lumi
        self.diffuse_flux_norm = diffuse_flux_norm
        self.diffuse_flux_ind = diffuse_flux_ind
        self.alerts, self.sources, self.uni_header = None, None, None
        self.timescale = kwargs.pop('timescale', 2.*86400.)
        self.sim_flux = None
        self.data_years = kwargs.pop('data_years', 1.)
        if self.diffuse_flux_norm is not None:
            self.N_per_dec_band() #initializes effective area
        self.manual_lumi = kwargs.pop('manual_lumi', 0.0)
        self.seed = kwargs.pop('seed', None)
        self.rng = np.random.RandomState(self.seed)

    def universe_firesong(self):
        return firesong_simulation('', density=self.density, Evolution=self.evolution,
                    Transient=True, timescale=self.timescale, fluxnorm = self.diffuse_flux_norm,
                    index=self.diffuse_flux_ind, LF = self.lumi, luminosity=self.manual_lumi, seed=self.seed)

    def create_universe(self):
        r'''
        Run FIRESONG for every year of data
        '''
        tmp_fls, tmp_dec, tmp_zs, tmp_tot = [],[],[],0.
        for i in range(int(self.data_years) / 1):
            uni = self.universe_firesong()
            tmp_dec.extend(uni['sources']['dec']), tmp_fls.extend(uni['sources']['flux'])
            tmp_zs.extend(uni['sources']['z'])
            tmp_tot += uni['total_flux']
        #Now do the fraction of a year
        if self.data_years % 1 != 0.0:
            uni = self.universe_firesong()
            add_src_num = int((self.data_years % 1) * len(uni['sources']['dec']))
            add_srcs_ind = self.rng.choice(list(range(len(uni['sources']['dec']))), add_src_num)
            tmp_dec.extend([uni['sources']['dec'][ind] for ind in add_srcs_ind])
            tmp_fls.extend([uni['sources']['flux'][ind] for ind in add_srcs_ind])
            tmp_zs.extend(uni['sources']['z'][ind] for ind in add_srcs_ind)
            tmp_tot += uni['total_flux'] * self.data_years % 1
        #fluxes are E^2 dN/dE at 100 TeV, convert now to dN/dE * DeltaT at 1 GeV
        tmp_fls = np.array(tmp_fls)
        #Time-integrated flux is over a duration of (1.+z) of the intrinsic burst time
        tmp_fls *= self.timescale * np.power(1e5, self.diffuse_flux_ind - 2.)*(1.+np.array(tmp_zs)) #*np.power(1, -1*self.diffuse_flux_ind)
        self.sources = {'dec': np.array(tmp_dec), 'flux': tmp_fls, 'z': np.array(tmp_zs)} 
        self.uni_header = uni['header']
        self.sim_flux = tmp_tot
        self.dec_band_from_decs()



def load_sig(cut = 'tight', stream = 'astro_numu'):
    r'''Load signalness distributions'''
    with open('/data/user/apizzuto/fast_response_skylab/alert_event_followup/signalness_distributions/{}_{}.csv'.format(stream, cut), 'r') as f:
        reader = csv.reader(f, delimiter=',')
        data = np.array(list(reader)).astype(float)
        sigs, heights = zip(*data)
    sigs, heights = list(sigs), list(heights)
    low_sig = find_nearest(np.linspace(0.025, 0.975, 20), sigs[0])
    sigs = np.linspace(low_sig, 0.975, len(heights))
    heights = np.maximum(heights, [0.0]*len(heights))
    return sigs, heights

def sample_from_hist(centers, heights, size = 1):
    r'''Given a histogram, sample from it, assuming
    uniform distributions within a bin'''
    cdf = np.cumsum(heights)
    cdf = cdf / cdf[-1]
    values = np.random.rand(size)
    value_bins = np.searchsorted(cdf, values)
    random_from_cdf = centers[value_bins]
    width = np.average(np.diff(centers))
    wiggle = np.random.uniform(low=- width / 2., high=width / 2., size=size)
    random_from_cdf += wiggle
    return random_from_cdf

def load_dec_dist(cut = 'gold'):
    if cut == 'gold':
        with open('/data/user/apizzuto/fast_response_skylab/alert_event_followup/signalness_distributions/Gold_backgrounds.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            data = np.array(list(reader)).astype(float)
            decs, heights = zip(*data)   
    else:
        with open('/data/user/apizzuto/fast_response_skylab/alert_event_followup/signalness_distributions/Gold_backgrounds.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            data = np.array(list(reader)).astype(float)
            decs, heights_gold = zip(*data) 
        with open('/data/user/apizzuto/fast_response_skylab/alert_event_followup/signalness_distributions/All_backgrounds.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            data = np.array(list(reader)).astype(float)
            decs, heights_all = zip(*data) 
        heights = np.array(heights_all) - np.array(heights_gold)
    decs = np.linspace(-0.9, 0.9, 10)
    heights = np.maximum(heights, [0.0]*len(heights))
    return decs, heights

def sample_declination(cut = 'gold', size = 1):
    sh = load_dec_dist(cut = cut)
    decs = sample_from_hist(sh[0], sh[1], size=size)
    return decs

def sample_signalness(stream='astro_numu', cut='tight', size = 1):
    r'''Load and sample from signalness distributions for various
    event streams'''
    if stream in ['astro_numu', 'conv_numu', 'conv_mu']:
        sh = load_sig(stream=stream, cut=cut)
    elif stream in ['signal', 'background', 'total']:
        if stream == 'signal':
            sh = load_sig(stream='astro_numu', cut = cut)
        elif stream == 'background':
            numu = load_sig(stream='conv_numu', cut=cut)
            mu = load_sig(stream='conv_mu', cut=cut)
            sh = [None, None]
            sh[0] = mu[0]
            sh[1] = numu[1] + mu[1]
        else:
            astro_numu = load_sig(stream='astro_numu', cut=cut)
            numu = load_sig(stream='conv_numu', cut=cut)
            mu = load_sig(stream='conv_mu', cut=cut)
            sh = [None, None]
            sh[0] = mu[0]
            sh[1] = numu[1] + mu[1] + astro_numu[1]
    else:
        print("Invalid stream name. Later, hater")
        return None
    sigs = sample_from_hist(sh[0], sh[1], size = size)
    return sigs

def load_alert_effA(stream = 'HESE', level='bronze'):
    r'''Load tabulated effective areas for different alert streams'''
    if stream == 'HESE':
        return load_HESE_effA(level = level)
    elif stream == 'EHE' or stream == 'GFU':
        return load_EHE_effA(level = level)
    else:
        print("Invalid alert stream")
        return None

def load_HESE_effA(level = 'bronze'):
    r'''Load starting events effective areas'''
    hese_tmp = np.load(eff_area_path + 'realtimeHESEv2_effA.npy').item()
    hese_ret = {}
    for key, val in hese_tmp.items():
        if level in key:
            ens = np.power(10., centers(val[1])) 
            effs = val[0] 
            hese_ret[key] = (ens, effs)
    return hese_ret

def load_EHE_effA(level = 'bronze'):
    r'''Load through-going events effective areas (this is actually
    GFU and EHE)'''
    ehe_tmp = np.load(eff_area_path + 'through_going_v2_alerts.npy').item()
    ehe_ret = {}
    for key, val in ehe_tmp.items():
        if level in key:
            ens = np.power(10., centers(val[1]))
            effs = val[0]
            ehe_ret[key] = (ens, effs)
    return ehe_ret

def aeff(en, dec, stream='HESE', level='bronze'):
    r'''Return effective area in m^2'''
    eff_A = load_alert_effA(stream=stream, level=level)
    for decs in eff_A.keys():
        d_range = ast.literal_eval(decs.split('_')[0])
        if dec > d_range[0] and dec <= d_range[1]:
            effA = eff_A[decs]
        else:
            pass
    aeffs = np.interp(en, effA[0], effA[1])
    aeffs = np.where(en > effA[0][0], aeffs, 0.)
    return aeffs

def spectrum(energy, gamma = -2., flux_norm = 1., ref_e = 1., cutoff = False):
    r'''Spectrum in units of GeV^-1 cm^-2 s^-1'''
    if cutoff and cutoff is not None:
        return flux_norm * np.power(energy/ref_e, gamma) \
             * np.exp(-1. * energy / cutoff)
    else:
        return flux_norm * np.power(energy / ref_e, gamma)
    
def calc_dnde(ens, dec, stream='HESE', level='bronze', gamma=-2.,
                 flux_norm=1., ref_e=1., cutoff=False):
    r'''Units of GeV^-1 s^-1'''
    eff = aeff(ens, dec, stream=stream, level=level)
    spec = spectrum(ens, gamma=gamma, flux_norm=flux_norm, ref_e=ref_e,
                       cutoff=cutoff)
    return eff*spec*1e4 #for m^2 to cm^2 in effective area

def calc_expected_signal_binned(ens, dec, stream='HESE', level='bronze', gamma=-2.,
                 flux_norm=1., ref_e=1., cutoff=False):
    r'''Units of s^-1'''
    eff = aeff(centers(ens), dec, stream=stream, level=level)
    spec = spectrum(centers(ens), gamma=gamma, flux_norm=flux_norm, ref_e=ref_e,
                       cutoff=cutoff)
    return eff*spec*np.diff(ens)*1e4

def calc_mean_n(ens, dec, stream='HESE', level='bronze', gamma=-2.,
                 flux_norm=1., ref_e=1., cutoff=False):
    return np.sum(calc_expected_signal_binned(ens, dec, stream=stream, level=level,
                          gamma=gamma,flux_norm=flux_norm, ref_e=ref_e,cutoff=cutoff))

def sample_signal_events(ens, dec, stream='HESE', level='bronze', gamma=-2.,
                 flux_norm=1., ref_e=1., cutoff=False):
    mu = calc_mean_n(ens, dec, stream=stream, level=level,
                          gamma=gamma,flux_norm=flux_norm, ref_e=ref_e,cutoff=cutoff)
    N = np.random.poisson(lam = mu)
    cut = 'tight' if level == 'gold' else 'loose'
    sigs = sample_signalness(stream='astro_numu', cut=cut, size=N)
    return N, sigs

def centers(arr):
    return arr[:-1] + np.diff(arr) / 2.

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_nearest_ind(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

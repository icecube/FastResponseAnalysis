import numpy as np
from glob import glob
import pandas as pd
import pickle
import csv
import ast
#from skylab.ps_injector import PointSourceInjector
#from skylab.datasets import Datasets
import sys
sys.path.append('/data/user/apizzuto/fast_response_skylab/alert_event_followup/FIRESONG/')
from Firesong import firesong_simulation

data_path = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/FIRESONG/Results/'
eff_area_path = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/effective_areas_alerts/'

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

bg_rates = {'HESE_gold': 0.4, 'HESE_bronze': 0.9, 'GFU_gold': 6.6, 'GFU_bronze': 13.8}

class TransientUniverse():
    r'''idk if bundling things like this will actually help out at all'''

    def __init__(self, lumi, evol, density, diffuse_flux_norm, diffuse_flux_ind,
                    **kwargs):

        self.evolution = evol
        self.density = density
        self.lumi = lumi
        self.diffuse_flux_norm = diffuse_flux_norm
        self.diffuse_flux_ind = diffuse_flux_ind
        self.alerts, self.sources, self.uni_header = None, None, None
        self.timescale = 2.*86400.
        self.sim_flux = None
        self.data_years = kwargs.pop('data_years', 1.)
        self.N_per_dec_band()

    def find_alerts(self):
        if self.uni_header is None:
            self.create_universe()
        signal_alerts = self.find_signal_alerts()
        background_alerts = self.find_background_alerts()
        alerts = signal_alerts + background_alerts
        self.alerts = alerts

    ##############################################3
    #### CHECK THE FLUXES, PROBABLY NEED TO DO SOME MULTIPLICATION
    ################################################
    def create_universe(self):
        r'''
        Run FIRESONG for every year of data
        '''
        tmp_fls, tmp_dec, tmp_tot = [],[],0.
        for i in range(int(self.data_years) / 1):
            uni = firesong_simulation('', density=self.density, Evolution=self.evolution,
                    Transient=True, timescale=self.timescale, fluxnorm = self.diffuse_flux_norm,
                    index=self.diffuse_flux_ind, LF = self.lumi)
            tmp_dec.extend(uni['sources']['dec']), tmp_fls.extend(uni['sources']['flux'])
            tmp_tot += uni['total_flux']
        #Now do the fraction of a year
        if self.data_years % 1 != 0.0:
            uni = firesong_simulation('', density=self.density, Evolution=self.evolution,
                    Transient=True, timescale=self.timescale, fluxnorm = self.diffuse_flux_norm,
                    index=self.diffuse_flux_ind, LF = self.lumi)
            add_src_num = int((self.data_years % 1) * len(uni['sources']['dec']))
            add_srcs_ind = np.random.choice(list(range(len(uni['sources']['dec']))), add_src_num)
            tmp_dec.extend([uni['sources']['dec'][ind] for ind in add_srcs_ind])
            tmp_fls.extend([uni['sources']['flux'][ind] for ind in add_srcs_ind])
            tmp_tot += uni['total_flux'] * self.data_years % 1
        #fluxes are E^2 dN/dE at 100 TeV, convert now to dN/dE * DeltaT at 1 GeV
        tmp_fls = np.array(tmp_fls)
        tmp_fls *= self.timescale * np.power(1, -1*self.diffuse_flux_ind)*np.power(1e5, self.diffuse_flux_ind - 2.)
        self.sources = {'dec': np.array(tmp_dec), 'flux': tmp_fls} 
        self.uni_header = uni['header']
        self.sim_flux = tmp_tot
        self.dec_band_from_decs()

    def find_background_alerts(self):
        bg_alerts = []
        for sample in ['HESE', 'GFU']:
                for cut, level in [('gold', 'tight'), ('bronze', 'loose')]:
                    N = np.random.poisson(lam = bg_rates[sample + '_' + cut])
                    sigs = sample_signalness(cut=level, stream='background', size=N)
                    bg_alerts.append((N, sigs, sample + '_' + cut, 'bg'))
        self.bg_alerts = bg_alerts
        return bg_alerts

    def dec_band_from_decs(self):
        bands = np.digitize(self.sources['dec'], bins = [-90., -5., 30., 90.]) - 1.
        self.sources['dec_bands'] = bands

    def N_per_dec_band(self):
        ens = np.logspace(1., 8., 501)
        tmp = {}
        for sample in ['HESE', 'GFU']:
            for level in ['gold', 'bronze']:
                tmp[sample + '_' + level] = []
                for dec in [-45., 0., 45]:
                    test_fl = 3.0 if sample == 'HESE' else 1.0 #eff A for HESE is all flavor
                    tmp[sample + '_' + level].append(calc_mean_n(ens, dec, stream=sample, 
                        level=level, gamma=self.diffuse_flux_ind*-1., flux_norm=test_fl))
        self.n_per_dec = tmp

    def find_signal_alerts(self):
        sig_alerts = [(0,0.,'','')]*len(self.sources['dec'])
        ens = np.logspace(1., 8., 501)
        for jjj, (src_dec, src_flux, src_bnd) in enumerate(zip(self.sources['dec'], self.sources['flux'], self.sources['dec_bands'])):
            for stream in ['GFU', 'HESE']:
                for cut, lev in [('gold', 'tight'), ('bronze', 'loose')]:
                    N = np.random.poisson(lam=self.n_per_dec[stream + '_' + cut][int(src_bnd)] * src_flux)
                    if N != 0.0:
                        sigs = sample_signalness(cut=lev, stream='signal', size = N)
                        sig_alerts[jjj] = (N, sigs, stream + '_' + cut,'sig')
        self.sig_alerts = sig_alerts
        return sig_alerts


# class Alert():

#     def __init__():
#         self.stream = 'HESE'
#         self.level = 'gold'
#         self.signal = False
#         self.dec = 0.0
#         self.skymap = None
#         self.time = None
#         self.additional_events = None

def load_sig(cut = 'tight', stream = 'astro_numu'):
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
    cdf = np.cumsum(heights)
    cdf = cdf / cdf[-1]
    values = np.random.rand(size)
    value_bins = np.searchsorted(cdf, values)
    random_from_cdf = centers[value_bins]
    wiggle = np.random.uniform(low=-0.025, high=0.025, size=size)
    random_from_cdf += wiggle
    return random_from_cdf


def sample_signalness(stream='astro_numu', cut='tight', size = 1):
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



# class BackgroundAlert(Alert):
#     def __init__(self):
#         print 'hi'

# class SignalAlert(Alert):
#     r'''TransientUniverse gets a list of Alerts maybe?'''

    # def __init__(self, dec, flux, timescale = 2.*86400.):
    #     self.dec, self.flux, self.timescale = dec, flux, timescale
    #     self.signalness = self.sample_signalness(stream='')




def load_alert_effA(stream = 'HESE', level='bronze'):
    if stream == 'HESE':
        return load_HESE_effA(level = level)
    elif stream == 'EHE' or stream == 'GFU':
        return load_EHE_effA(level = level)
    else:
        print("Invalid alert stream")
        return None

def load_HESE_effA(level = 'bronze'):
    hese_tmp = np.load(eff_area_path + 'realtimeHESEv2_effA.npy').item()
    hese_ret = {}
    for key, val in hese_tmp.items():
        if level in key:
            ens = np.power(10., centers(val[1])) 
            effs = val[0] 
            hese_ret[key] = (ens, effs)
    return hese_ret

def load_EHE_effA(level = 'bronze'):
    ehe_tmp = np.load(eff_area_path + 'through_going_v2_alerts.npy').item()
    ehe_ret = {}
    for key, val in ehe_tmp.items():
        if level in key:
            ens = np.power(10., centers(val[1]))
            effs = val[0]
            ehe_ret[key] = (ens, effs)
    #COMMENTED OUT BLOCK FROM BEFORE I FOUND MC FILE
    #decs = ['[-90.0, -5.0]','[-5.0, 30.0]','[30.0, 90.0]']
    #for dec in decs:
    #    dec_list = ast.literal_eval(dec)
    #    with open(eff_area_path + 'EHE_GFU_{}_{}.csv'.format(int(dec_list[0]), int(dec_list[1])), 'r') as f:
    #        reader = csv.reader(f, delimiter=',')
    #        data = np.array(list(reader)).astype(float)
    #        ens, effs = zip(*data)
    #    ehe_ret[dec+'_'+level]=(np.array(ens)*1e3, np.array(effs))
    #return ehe_ret
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


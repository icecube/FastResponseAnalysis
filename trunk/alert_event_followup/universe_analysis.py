import numpy as np
from glob import glob
import pandas as pd
import scipy.stats as st
import pickle
import csv
import ast
import sys
sys.path.append('/data/user/apizzuto/fast_response_skylab/alert_event_followup/FIRESONG/')
from Firesong import firesong_simulation
from transient_universe import TransientUniverse, SteadyUniverse

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

bg_trials = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/bg/'
signal_trials = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/fits/'

class UniverseAnalysis():
    r'''Given cosmological parameters, calculate the expected TS distribution
        from triggering short timescale analyses on alert events'''

    def __init__(self, lumi, evol, density, diffuse_flux_norm, diffuse_flux_ind,
                    **kwargs):
        self.lumi = lumi
        self.evol = evol
        self.density = density
        self.diffuse_flux_norm = diffuse_flux_norm
        self.diffuse_flux_ind = diffuse_flux_ind
        self.deltaT = kwargs.pop('deltaT', None)
        self.transient = True if self.deltaT is not None else False
        self.seed = kwargs.pop('seed', 1234)
        if self.transient:
            self.universe = TransientUniverse(self.lumi, self.evol, self.density,
                self.diffuse_flux_norm, self.diffuse_flux_ind, seed=self.seed, **kwargs)
        else:
            self.universe = SteadyUniverse(self.lumi, self.evol, self.density,
                self.diffuse_flux_norm, self.diffuse_flux_ind, seed=self.seed, **kwargs)
        self.smear = kwargs.pop('smeared', True)
        self.smear_str = 'smeared/' if self.smear else 'norm_prob/'
        self.verbose = kwargs.pop('verbose', False)
        self.initialize_universe() 

    def initialize_universe(self):
        if self.verbose:
            print("Simulating universe with specified cosmological parameters")
        self.universe.create_universe()
        self.universe.find_alerts()
        self.universe.find_alert_skymaps()
        self.universe.additional_signal_events()

    def make_alerts_dataframe(self):
        alerts = {'signalness': [], 'declination': [], 'background': [], 
          'skymap_ind': [], 'stream': [], 'skymap_dec': [],
         'extra_evs': []}
        for k in self.universe.bg_alerts.keys():
            if self.universe.bg_alerts[k][0] > 0:
                alerts['signalness'].extend(self.universe.bg_alerts[k][1])
                alerts['declination'].extend(self.universe.bg_alerts[k][2])
                alerts['background'].extend([True]*self.universe.bg_alerts[k][0])
                alerts['skymap_ind'].extend(self.universe.bg_alerts[k][4]) 
                alerts['skymap_dec'].extend(self.universe.bg_alerts[k][3]) 
                alerts['stream'].extend([k]*self.universe.bg_alerts[k][0])
                alerts['extra_evs'].extend([0]*self.universe.bg_alerts[k][0])
        for k in self.universe.sig_alerts.keys():
            for jj in range(len(self.universe.sig_alerts[k])):
                if self.universe.sig_alerts[k][jj][0] == 0:
                    continue
                else:
                    alerts['signalness'].append(self.universe.sig_alerts[k][jj][1][0])
                    alerts['declination'].append(np.radians(self.universe.sources['dec'][jj]))
                    alerts['background'].append(False)
                    alerts['skymap_ind'].append(self.universe.skymaps[k][jj][1])
                    alerts['skymap_dec'].append(self.universe.skymaps[k][jj][0])
                    alerts['stream'].append(k)
                    alerts['extra_evs'].append(self.universe.extra_events[k][jj])
        alerts = pd.DataFrame(alerts)
        self.alert_df = alerts

    def reinitialize_universe(self):
        if self.verbose:
            print("Recreating universe for more trials, updating seed")
        self.seed = 1 if self.seed is None else self.seed + 1
        self.universe.seed = self.seed
        self.universe.create_universe()
        self.universe.find_alerts()
        self.universe.find_alert_skymaps()
        self.universe.additional_signal_events()

    def calculate_ts(self, only_gold = False, calc_p=True):
        ts, sigs, ps = [], [], []
        self.alert_df['TS'] = [None] * len(self.alert_df['background'])
        self.alert_df['pval'] = [None] * len(self.alert_df['background'])
        for index, alert in self.alert_df.iterrows():
            if alert['background']:
                if calc_p:
                    t, p = self.background_alert_trials(alert['skymap_ind'], calc_p=calc_p)
                    ts.append(t); ps.append(p)
                else:
                    ts.append(self.background_alert_trials(alert['skymap_ind'], calc_p=calc_p))
                sigs.append(alert['signalness'])
                self.alert_df.loc[self.alert_df.index == index, 'TS'] = ts[-1]
                if calc_p:
                    self.alert_df.loc[self.alert_df.index == index, 'pval'] = ps[-1]
            else:
                if calc_p:
                    t, p = self.signal_alert_trials(alert['skymap_ind'], alert['extra_evs'], calc_p=calc_p)
                    ts.append(t); ps.append(p)
                else:
                    ts.append(self.signal_alert_trials(alert['skymap_ind'], alert['extra_evs'], calc_p=calc_p))
                sigs.append(alert['signalness'])
                self.alert_df.loc[self.alert_df.index == index, 'TS'] = ts[-1]
                if calc_p:
                    self.alert_df.loc[self.alert_df.index == index, 'pval'] = ps[-1]
        ts, sigs = np.array(ts), np.array(sigs)
        if only_gold:
            gold = []
            for index, alert in self.alert_df.iterrows():
                if 'gold' in alert['stream']:
                    gold.append(True)
                else:
                    gold.append(False)
            gold = np.array(gold)
            ts, sigs = ts[gold], sigs[gold]
        TS = np.sum(sigs * ts) / sigs.size
        self.TS = TS
        return TS
        
    def background_alert_trials(self, ind, calc_p=True):
        if self.transient:
            trials_file = glob(bg_trials + self.smear_str + 'index_{}_*_time_{:.1f}.pkl'.format(ind, self.deltaT))[0]
            trials = np.load(trials_file)
            ts = np.random.choice(trials['ts_prior'])
            if calc_p:
                if ts == 0:
                    pval = 1.0
                else:
                    pval = float(np.count_nonzero(np.array(trials['ts_prior']) >= ts)) / np.array(trials['ts_prior']).size
                    if pval == 0.:
                        pval = 1./np.array(trials['ts_prior']).size
        else:
            fs = glob(bg_trials + self.smear_str + 'index_{}_steady_seed_*.pkl'.format(ind))
            #f = np.random.choice(fs)
            #trials = np.load(f) 
            tmp_ds = [np.load(f) for f in fs]
            trials = {k:[] for k in tmp_ds[0].keys()}
            for k in trials.keys():
                for d in tmp_ds:
                    trials[k] += d[k] 
            for k in trials.keys():
                trials[k] = np.array(trials[k])
            ts = np.random.choice(trials['TS'])
            if calc_p:
                if ts == 0:
                    pval = 1.0
                else:
                    pval = float(np.count_nonzero(trials['TS'] >= ts)) / trials['TS'].size
                    if pval == 0.:
                        pval = 1./np.array(trials['ts_prior']).size
        #ts = np.random.choice(trials['ts_prior'])
        del trials
        if calc_p:
            return ts, pval
        else:
            return ts

    def signal_alert_trials(self, ind, N, calc_p = True):
        if N == 0:
            ts = self.background_alert_trials(ind, calc_p = False)
        else:
            if self.transient:
                trials_file = glob(signal_trials + self.smear_str + 'index_{}_*_time_{:.1f}.pkl'.format(ind, self.deltaT))[0]
                trials = np.load(trials_file)
            else:
                fs = glob(signal_trials + self.smear_str + 'index_{}_steady_seed_*.pkl'.format(ind))
                t_file = np.random.choice(fs)
                trials = np.load(t_file)
                #trials = np.load(signal_trials + 'index_{}_steady.pkl'.format(ind))
            ns_key = 'true_ns' if self.transient else 'inj_nsignal'
            if N <= 10:
                inds = np.argwhere(np.array(trials[ns_key]) == N).flatten()
                if len(inds) == 0:
                    inds = np.argwhere(np.abs(np.array(trials[ns_key]) - N) < 4).flatten()
                    if len(inds) == 0:
                        if self.verbose:
                            print("No trials near {}".format(N))
                        inds = np.argmin(np.abs(np.array(trials[ns_key]) - N)).flatten()    
            else:
                inds = np.argwhere(np.abs(np.array(trials[ns_key]) - N) < 10).flatten()
                if len(inds) == 0:
                    if self.verbose:
                        print("NO TRIALS WITH {} INJECTED EVENTS".format(N))
                    inds = np.argwhere(np.array(trials[ns_key]) == np.max(trials[ns_key])).flatten()
            ts = np.array(trials['ts'])[inds] if self.transient else np.array(trials['TS'])[inds]
            ts = np.random.choice(ts)
            del trials
        if calc_p:
            pval = self.calculate_trial_pvalue(ind, ts)
            return ts, pval
        else:
            return ts

    def calculate_trial_pvalue(self, ind, TS):
        if TS == 0:
            return 1.
        if self.transient:
            trials_file = glob(bg_trials + self.smear_str + 'index_{}_*_time_{:.1f}.pkl'.format(ind, self.deltaT))[0]
            trials = np.load(trials_file)
            pval = float(np.count_nonzero(np.array(trials['ts_prior']) >= TS)) / np.array(trials['ts_prior']).size
            if pval == 0.:
                pval = 1./np.array(trials['ts_prior']).size
        else:
            fs = glob(bg_trials + self.smear_str + 'index_{}_steady_seed_*.pkl'.format(ind))
            tmp_ds = [np.load(f) for f in fs]
            trials = {k:[] for k in tmp_ds[0].keys()}
            for k in trials.keys():
                for d in tmp_ds:
                    trials[k] += d[k]
            for k in trials.keys():
                trials[k] = np.array(trials[k])
            pval = float(np.count_nonzero(trials['TS'] >= TS)) / trials['TS'].size
            if pval == 0.:
                pval = 1./trials['TS'].size
            del trials
        return pval

    def calculate_binomial_pvalue(self, only_gold=False):
        if self.TS is None:
            self.calculate_ts(only_gold = only_gold, calc_p=True) 
        plist = self.alert_df['pval']
        if only_gold:
            stream_msk = self.alert_df['stream']
            stream_msk = ['gold' in al for al in stream_msk] 
            plist = plist[stream_msk]
        obs_p = 1.
        plist = sorted(plist)
        for i, p in enumerate(plist):
            
            tmp = st.binom_test(i+1, len(plist), p, alternative='greater')
            if tmp < obs_p and tmp != 0.0:
                if tmp == 0.0:
                    print("WHY DOES THE BINOMIAL VALUE EQUAL ZERO")
                obs_p = tmp
        self.binom_p = obs_p
        return obs_p 

import numpy as np
from glob import glob
import pandas as pd
import pickle
import csv
import ast
import sys
sys.path.append('/data/user/apizzuto/fast_response_skylab/alert_event_followup/FIRESONG/')
from Firesong import firesong_simulation
from transient_universe import TransientUniverse

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
        self.seed = kwargs.pop('seed', None)
        if self.transient:
            self.universe = TransientUniverse(self.lumi, self.evol, self.density,
                self.diffuse_flux_norm, self.diffuse_flux_ind, seed=self.seed, **kwargs)
        else:
            self.universe = SteadyUniverse(self.lumi, self.evol, self.density,
                self.diffuse_flux_norm, self.diffuse_flux_ind, seed=self.seed, **kwargs)
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
                if self.universe.sig_alerts[k][jj] == (0, 0.0):
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

    def calculate_ts(self):
        ts, sigs = [], []
        self.alert_df['TS'] = [None] * len(self.alert_df['background'])
        for index, alert in self.alert_df.iterrows():
            if alert['background']:
                ts.append(self.background_alert_trials(alert['skymap_ind']))
                sigs.append(alert['signalness'])
                self.alert_df.loc[self.alert_df.index == index, 'TS'] = ts[-1]
            else:
                ts.append(self.signal_alert_trials(alert['skymap_ind'], alert['extra_evs']))
                sigs.append(alert['signalness'])
                self.alert_df.loc[self.alert_df.index == index, 'TS'] = ts[-1]
        ts, sigs = np.array(ts), np.array(sigs)
        TS = np.sum(sigs * ts) / sigs.size
        self.TS = TS
        return TS
        
    def background_alert_trials(self, ind):
        if ind == 88:
            return np.nan
        if self.transient:
            trials = np.load(bg_trials + 'index_{}_time_{:.1f}.pkl'.format(ind, self.deltaT))
        else:
            trials = np.load(bg_trials + 'index_{}_steady.pkl'.format(ind))
        ts = np.random.choice(trials['ts_prior'])
        del trials
        return ts

    def signal_alert_trials(self, ind, N):
        if ind == 88:
            return np.nan
        if N == 0:
            ts = self.background_alert_trials(ind)
        else:
            if self.transient:
                trials = np.load(signal_trials + 'index_{}_time_{:.1f}.pkl'.format(ind, self.deltaT))
            else:
                trials = np.load(signal_trials + 'index_{}_steady.pkl'.format(ind))
            if N <= 10:
                inds = np.argwhere(np.array(trials['true_ns']) == N).flatten()
                if len(inds) == 0:
                    inds = np.argwhere(np.abs(np.array(trials['true_ns']) - N) < 4).flatten()
                    if len(inds) == 0:
                        if self.verbose:
                            print("No trials near {}".format(N))
                        inds = np.argmin(np.abs(np.array(trials['true_ns']) - N)).flatten()    
            else:
                inds = np.argwhere(np.abs(np.array(trials['true_ns']) - N) < 10).flatten()
                if len(inds) == 0:
                    if self.verbose:
                        print("NO TRIALS WITH {} INJECTED EVENTS".format(N))
                    inds = np.argwhere(np.array(trials['true_ns']) == np.max(trials['true_ns'])).flatten()
            ts = np.array(trials['ts'])[inds]
            ts = np.random.choice(ts)
            del trials
        return ts
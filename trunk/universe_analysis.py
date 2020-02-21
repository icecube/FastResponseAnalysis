import numpy as np
from glob import glob
import pandas as pd
import pickle
import csv
import ast
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

bg_trials = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/bg/'
signal_trials = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/fits/'

def UniverseAnalysis():
    r'''Given cosmological parameters, calculate the expected TS distribution
        from triggering short timescale analyses on alert events'''

    def __init__(self, ):
        self.something = None
        self.deltaT = 2. * 86400.
        self.universe = TransientUniverse(lumi, evol, density,
            diffuse_flux_norm, diffuse_flux_ind, **kwargs)
        self.verbose = kwargs.pop('verbose', False)

    def initialize_universe(self):
        if self.verbose:
            print("Simulating universe with specified cosmological parameters")
        self.universe.create_universe()
        self.universe.find_alerts()

    def background_trials(self, N_trials=100):
        return 0

    def signal_trials(self, N_trials=100):
        return 0.




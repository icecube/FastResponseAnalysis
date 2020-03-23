#!/usr/bin/env python

import os, sys, time, pickle
sys.path.append('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/')
import numpy as np
from transient_universe import TransientUniverse
from universe_analysis import UniverseAnalysis
import argparse

parser = argparse.ArgumentParser(description='Calculate TS distributions')
parser.add_argument('--n', type=int, default=1000,
                        help = 'Number of trials')
parser.add_argument('--density', type=float, default=1e-9,
                        help = 'Local source density')
parser.add_argument('--LF', type=str, default='SC', help ='luminosity function')
parser.add_argument('--evol', type=str, default='HB2006SFR', help='Evolution')
#parser.add_argument('--gold', default=False, action='store_true', help='Only analyze gold alerts')
args = parser.parse_args()

TS = []
TS_gold = []

density = args.density
evol = args.evol
lumi = args.LF
#only_gold = args.gold

uni = UniverseAnalysis(lumi, evol, density, 1.01e-8, 2.19, deltaT=2*86400., 
        data_years=8)
uni.initialize_universe()
uni.make_alerts_dataframe()
TS.append(uni.calculate_ts())

for jj in range(args.n - 1):
    uni.reinitialize_universe()
    uni.make_alerts_dataframe()
    TS.append(uni.calculate_ts(only_gold = False))
    TS_gold.append(uni.calculate_ts(only_gold = True))
        
#TS = np.array(TS)
#TS_gold = np.array(TS_gold)
TS = np.array([TS, TS_gold])

np.save('/data/user/apizzuto/fast_response_skylab/alert_event_followup/ts_distributions/ts_dists_5year_density_{:.2e}_evol_{}_lumi_{}.npy'.format(density, evol, lumi), TS)

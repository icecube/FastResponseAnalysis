#!/usr/bin/env python

import numpy as np
import os, sys, argparse, pickle
from astropy.time import Time
from astropy.time import TimeDelta
from numpy.lib.recfunctions import append_fields

from fast_response.FastResponseAnalysis import FastResponseAnalysis

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--sinDec', type=float,default=None,
                    help='sinDec of source')
parser.add_argument('--deltaT', type=float, default=None,
                    help='Time Window in seconds')
parser.add_argument('--month', type=int, default=6,
                    help='Month for the analysis in 2018 - 2019, for seasonal variation stuff')
parser.add_argument('--ntrials', type=int, default=100000,
                    help='Number of trials to run')
parser.add_argument('--rng_seed', type=int, default=1,
                    help='Random number seed for repeated jobs')
args = parser.parse_args()

deltaT = args.deltaT / 86400.

if args.month < 1 or args.month > 12:
    print("Please enter a valid month (1-12). Exiting")
    sys.exit()
if args.month <= 6:
    start = "2017-{:02d}-01 00:00:00".format(args.month)
else:
    start = "2016-{:02d}-01 00:00:00".format(args.month)

stop = (Time(start) + TimeDelta(deltaT)).iso
dec = np.arcsin(args.sinDec)

f = FastResponseAnalysis("0., {}".format(dec*180. / np.pi), start, stop, Name="test", save=False, seed=args.rng_seed)

results = f.llh.do_trials(args.ntrials, src_ra = 0., src_dec = dec)
#names = results.dtype.names
#names = list(names)
#names.remove('spectrum')
#results = results[names]

bg_trials = np.array(results['TS'])
bg_trials = {'n_zero': np.count_nonzero(bg_trials == 0.),
                'tsd': bg_trials[bg_trials != 0.0]}
with open('/data/user/apizzuto/fast_response_skylab/fast-response/fast_response/analysis_checks/bg/ps_sinDec_{}_deltaT_{}_month_{}_seed_{}.pkl'.format(args.sinDec, args.deltaT, args.month, args.rng_seed), 'w') as f:
    pickle.dump(bg_trials, f)


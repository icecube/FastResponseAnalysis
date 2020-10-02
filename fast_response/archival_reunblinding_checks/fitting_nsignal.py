#!/usr/bin/env python

import numpy as np
import os, sys, argparse
from astropy.time import Time
from astropy.time import TimeDelta
from numpy.lib.recfunctions import append_fields

from fast_response.FastResponseAnalysis import FastResponseAnalysis

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--sinDec', type=float,default=None,
                    help='sinDec of source')
parser.add_argument('--deltaT', type=float, default=None,
                    help='Time Window in seconds')
args = parser.parse_args()

gammas = np.linspace(2., 3., 3)
nsigs = np.linspace(1., 15., 15)
deltaT = args.deltaT / 86400.

start = "2017-06-01 00:00:00"
stop = (Time(start) + TimeDelta(deltaT)).iso
dec = np.arcsin(args.sinDec)

f = FastResponseAnalysis("0., {}".format(dec*180. / np.pi), start, stop, save=False)
results = None

for gamma in gammas:
    f = FastResponseAnalysis("0., {}".format(dec*180. / np.pi), start, stop, save=False)
    inj = f.initialize_injector(gamma=gamma)
    for nsig in nsigs:
        result = f.llh.do_trials(100, src_ra = 0., src_dec = dec, injector = inj, mean_signal=int(nsig), poisson=False)
        result = append_fields(result, 'gamma', [gamma]*len(result), usemask=False)
        if results == None:
            results = result
        else:
            results = np.append(results, result)

print(results)
np.save('/data/user/apizzuto/fast_response_skylab/fast-response/fast_response/analysis_checks/fit/nsignal_sinDec_{}_deltaT_{}.npy'.format(args.sinDec, args.deltaT),
            results)

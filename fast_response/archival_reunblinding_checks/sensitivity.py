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
parser.add_argument('--month', type=int, default=6,
                    help='Month for the analysis in 2018 - 2019, for seasonal variation stuff')
args = parser.parse_args()

gammas = np.linspace(2., 3., 3)
nsigs = np.linspace(1., 15., 15)
deltaT = args.deltaT / 86400.

if args.month < 1 or args.month > 12:
    print "Please enter a valid month (1-12). Exiting"
    sys.exit()
if args.month <= 6:
    start = "2017-{:02d}-01 00:00:00".format(args.month)
else:
    start = "2016-{:02d}-01 00:00:00".format(args.month)

stop = (Time(start) + TimeDelta(deltaT)).iso
dec = np.arcsin(args.sinDec)

f = FastResponseAnalysis("0., {}".format(dec*180. / np.pi), start, stop, save=False)
results = None

for gamma in gammas:
    f = FastResponseAnalysis("0., {}".format(dec*180. / np.pi), start, stop, save=False)
    inj = f.initialize_injector(gamma=gamma)
    for nsig in nsigs:
        result = f.llh.do_trials(1000, src_ra = 0., src_dec = dec, injector = inj, mean_signal=int(nsig), poisson=True)
        result = append_fields(result, 'gamma', [gamma]*len(result), usemask=False)
        result = append_fields(result, 'mean_ninj', [nsig]*len(result), usemask=False)
        result = append_fields(result, 'flux', [inj.mu2flux(int(nsig))]*len(result), usemask=False)
        names = result.dtype.names
        names = list(names)
        names.remove('spectrum')
        result = result[names]
        if results == None:
            results = result
        else:
            results = np.append(results, result)

print(results)
np.save('/data/user/apizzuto/fast_response_skylab/fast-response/fast_response/analysis_checks/sensitivity/nsignal_sinDec_{}_deltaT_{}_month_{}.npy'.format(args.sinDec, args.deltaT, args.month), results)

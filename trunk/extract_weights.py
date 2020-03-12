#!/usr/bin/env python

r'''Major discrepancy between old and new
    fast response is from the spatial and 
    energy weights. This script is just to 
    calculate a whole bunch of those and make
    scatter plots 

    Author: Alex Pizzuto
    Date: October 21, 2019
    '''

import os.path,sys,argparse
import subprocess
import warnings
warnings.filterwarnings("ignore")

base_path = os.path.join('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/','')
sys.path.append(base_path)

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--dec', type=float, default=None,
                    help='Dec in degrees or path to skymap')
parser.add_argument('--extension', type=float, default=None,
                    help="Source extension in degrees")
parser.add_argument('--deltaT', type=float, default=None,
                    help="Analysis duration in days (default=1)")
args = parser.parse_args()

import pickle
import logging as log
from astropy.time import Time
from astropy.coordinates import Angle
import astropy.units as u
from FastResponseAnalysis import FastResponseAnalysis

log.basicConfig(level=log.ERROR)
source = {}
source['Name'] = 'Weighting'
source['save'] = False
source['Extension'] = args.extension
source['Location'] = "0.0, {}".format(args.dec)
start = "2017-01-01 00:00:00"
mjd = 58119.0
duration = args.deltaT
results = []

for i in range(300):
    if duration > 5. and i > 100:
        continue
    start_mjd = mjd + (i*duration)
    stop_mjd = start_mjd + duration
    if (start_mjd < 58309.74448310) and (stop_mjd > 58309.74448310):
        continue
    tstart_str = Time(start_mjd, format='mjd', scale = 'utc').iso
    tstop_str = Time(stop_mjd, format='mjd', scale='utc').iso
    f = FastResponseAnalysis(source['Location'], tstart_str, tstop_str, **source)
    ts = f.unblind_TS()
    try:
        ns_profile = f.ns_scan()
    except:
        ns_profile = None
    p = f.calc_pvalue(ntrials = 10000)
    evs = f.coincident_events
    tsd = f.tsd
    results.append((ts, p, evs, tsd, ns_profile))
    

with open('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/analysis_checks/weights/weights_dec_{:.2f}_time_{:.2f}_ext_{}_results.pickle'.format(args.dec, args.deltaT, args.extension), 'wb') as f:
    pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)



#!/usr/bin/env python

import numpy as np
import healpy as hp
import os, sys, argparse, pickle
from astropy.time import Time
from astropy.time import TimeDelta
from numpy.lib.recfunctions import append_fields
from astropy.io import fits
from glob import glob

base_path = os.path.join('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/','')
sys.path.append(base_path)

from FastResponseAnalysis import FastResponseAnalysis

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--index', type=int,default=None,
                    help='skymap index')
parser.add_argument('--deltaT', type=float, default=None,
                    help='Time Window in seconds')
parser.add_argument('--ntrials', type=int, default = 10000,
                        help='Trials')
args = parser.parse_args()

output_paths = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/results/'

skymap_files = glob('/data/ana/realtime/alert_catalog_v2/2yr_prelim/fits_files/Run13*.fits.gz')
skymap_fits, skymap_header = hp.read_map(skymap_files[args.index], h=True, verbose=False)
skymap_header = {name: val for name, val in skymap_header}
ev_mjd = skymap_header['EVENTMJD']
ev_run, ev_id = skymap_header['RUNID'], skymap_header['EVENTID']
source = {"Skipped Events": [(ev_run, ev_id)]}

deltaT = args.deltaT / 86400.
event_mjd = ev_mjd
start_mjd = event_mjd - (deltaT / 2.)
stop_mjd = event_mjd + (deltaT / 2.)
start = Time(start_mjd, format='mjd').iso
stop = Time(stop_mjd, format='mjd').iso

f = FastResponseAnalysis(skymap_files[args.index], start, stop, save=False, 
                            alert_event=True, **source)
inj = f.initialize_injector(gamma=2.5) #just put this here to initialize f.spatial_prior
ts = f.unblind_TS()
res = f.save_results(alt_path = output_paths)


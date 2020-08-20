#!/usr/bin/env python

import numpy as np
import healpy as hp
import os, sys, argparse
from astropy.time import Time
from astropy.time import TimeDelta
from numpy.lib.recfunctions import append_fields
from astropy.io import fits
from glob import glob
import pickle

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
parser.add_argument('--smear', default=False, action='store_true',
                    help='Include systematics by smearing norm. prob.')
args = parser.parse_args()


#skymaps_path = '/data/user/steinrob/millipede_scan_archive/fits_v3_prob_map/'
#files = glob(skymaps_path + '*.fits')
output_paths = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/bg/'

skymap_files = glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run1*.fits.gz')

skymap_fits, skymap_header = hp.read_map(skymap_files[args.index], h=True, verbose=False)
skymap_header = {name: val for name, val in skymap_header}
ev_mjd = skymap_header['EVENTMJD']

run_id = skymap_header['RUNID']
event_id = skymap_header['EVENTID']

deltaT = args.deltaT / 86400.

#skymap_fits = fits.open(files[args.index])[0]
#event_mjd = skymap_fits.header['TIME_MJD']
event_mjd = ev_mjd    #58000.000 #HARDCODE SO THAT THERE IS REAL DATA
start_mjd = event_mjd - (deltaT / 2.)
stop_mjd = event_mjd + (deltaT / 2.)

start = Time(start_mjd, format='mjd').iso
stop = Time(stop_mjd, format='mjd').iso

trials_per_sig = args.ntrials

tsList_prior  = []
tsList        = []
nsList        = []
nsList_prior  = []
true_ns       = []
ra            = []
dec           = []

seed_counter = 0

f = FastResponseAnalysis(skymap_files[args.index], start, stop, save=False, 
                            alert_event=True, smear=args.smear, alert_type='track')
inj = f.initialize_injector(gamma=2.5) #just put this here to initialize f.spatial_prior
for jj in range(trials_per_sig):
    seed_counter += 1
    try:
        val = f.llh.scan(0.0, 0.0, scramble=True, seed = seed_counter,
                spatial_prior = f.spatial_prior, time_mask = [deltaT / 2., event_mjd],
                pixel_scan = [f.nside, 4.0], inject = None)
        tsList_prior.append(val['TS_spatial_prior_0'].max())
        tsList.append(val['TS'].max())
        max_prior   = np.argmax(val['TS_spatial_prior_0'])
        max_noPrior = np.argmax(val['TS'])
        nsList_prior.append(val['nsignal'][max_prior])
        nsList.append(val['nsignal'][max_noPrior])
        ra.append(val['ra'][max_prior])
        dec.append(val['dec'][max_prior])
    except ValueError:
        tsList_prior.append(0.0)
        tsList.append(0.0)
        nsList_prior.append(0.0)
        nsList.append(0.0)
        ra.append(0.0)
        dec.append(0.0)

results = {'ts_prior': tsList_prior, 'ts': tsList, 'ns_prior': nsList_prior,
            'ns': nsList, 'ra': ra, 'dec': dec}

smear_str = 'smeared/' if args.smear else 'norm_prob/'
with open(output_paths + smear_str + 'index_{}_run_{}_event_{}_time_{}.pkl'.format(args.index, run_id, event_id, args.deltaT), 'w') as fi:
    pickle.dump(results, fi, protocol=pickle.HIGHEST_PROTOCOL)

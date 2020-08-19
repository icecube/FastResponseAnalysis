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
parser.add_argument('--ntrials', type=int, default = 100,
                        help='Trials per injection strength')
parser.add_argument('--smear', default=False, action='store_true',
                    help='Include systematics by smearing norm. prob.')
args = parser.parse_args()

#skymaps_path = '/data/user/steinrob/millipede_scan_archive/fits_v3_prob_map/'
#files = glob(skymaps_path + '*.fits')
output_paths = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/fits/'

skymap_files = glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run1*.fits.gz')

skymap_fits, skymap_header = hp.read_map(skymap_files[args.index], h=True, verbose=False)
skymap_header = {name: val for name, val in skymap_header}
ev_mjd = skymap_header['EVENTMJD']
run_id = skymap_header['RUNID']
event_id = skymap_header['EVENTID']

gammas = [2.5] #np.linspace(2., 3., 3)
nsigs = [1., 2., 3., 4., 5, 6., 7., 8., 9., 10., 15., 20., 25., 30., 35., 40.]
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
gammaList     = []
mean_ninj     = []
flux_list     = []
true_ras      = []
true_decs     = []

seed_counter = 0

for gamma in gammas:
    f = FastResponseAnalysis(skymap_files[args.index], start, stop, save=False, 
                alert_event=True, smear=args.smear)
    inj = f.initialize_injector(gamma=gamma)
    scale_arr = []
    step_size = 10
    for i in range(1,20*step_size + 1, step_size):
        scale_arr.append([])
        for j in range(5):
            scale_arr[-1].append(inj.sample(i, poisson=False)[0][0])
    scale_arr = np.median(scale_arr, axis=1)
    try:
        scale_factor = np.min(np.argwhere(scale_arr > 0))*step_size + 1.
    except:
        print("Scale factor thing for prior injector didn't work")
        scale_factor = 1.
    for nsig in nsigs:
        for jj in range(trials_per_sig):
            seed_counter += 1
            ni, sample, true_ra, true_dec = inj.sample(nsig*scale_factor, poisson = False, return_position = True)
            if sample is not None:
                sample['time'] = event_mjd
            try:
                val = f.llh.scan(0.0, 0.0, scramble=True, seed = seed_counter,
                        spatial_prior = f.spatial_prior, time_mask = [deltaT / 2., event_mjd],
                        pixel_scan = [f.nside, 4.0], inject = sample)
                tsList_prior.append(val['TS_spatial_prior_0'].max())
                tsList.append(val['TS'].max())
                max_prior   = np.argmax(val['TS_spatial_prior_0'])
                max_noPrior = np.argmax(val['TS'])
                nsList_prior.append(val['nsignal'][max_prior])
                nsList.append(val['nsignal'][max_noPrior])
                true_ns.append(val['n_inj'][max_prior])
                ra.append(val['ra'][max_prior])
                dec.append(val['dec'][max_prior])
                gammaList.append(gamma)
                mean_ninj.append(nsig)
                flux_list.append(inj.mu2flux(nsig))
                true_ras.append(true_ra[0])
                true_decs.append(true_dec[0])
            except ValueError:
                tsList_prior.append(0.0)
                tsList.append(0.0)
                nsList_prior.append(0.0)
                nsList.append(0.0)
                true_ns.append(0.0)
                ra.append(0.0)
                dec.append(0.0)
                gammaList.append(gamma)
                mean_ninj.append(nsig)
                flux_list.append(inj.mu2flux(nsig*scale_factor))
                true_ras.append(true_ra[0])
                true_decs.append(true_dec[0])

results = {'ts_prior': tsList_prior, 'ts': tsList, 'ns_prior': nsList_prior,
            'ns': nsList, 'true_ns': true_ns, 'ra': ra, 'dec': dec,
            'gamma': gammaList, 'mean_ninj': mean_ninj, 'flux': flux_list,
            'true_ra': true_ras, 'true_dec': true_decs}

smear_str = 'smeared/' if args.smear else 'norm_prob/'
with open(output_paths + smear_str + 'index_{}_run_{}_event_{}_time_{}.pkl'.format(args.index, run_id, event_id, args.deltaT), 'w') as fi:
    pickle.dump(results, fi, protocol=pickle.HIGHEST_PROTOCOL)

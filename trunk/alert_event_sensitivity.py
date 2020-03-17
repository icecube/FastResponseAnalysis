#!/usr/bin/env python

import numpy as np
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
args = parser.parse_args()

skymaps_path = '/data/user/steinrob/millipede_scan_archive/fits_v3_prob_map/'
files = glob(skymaps_path + '*.fits')
output_paths = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/sensitivity_ts_distributions/'

gammas = [2.5] #np.linspace(2., 3., 3)
nsigs = [1., 2., 3., 4., 6., 8., 10., 15., 20., 25., 30., 50.]
deltaT = args.deltaT / 86400.

skymap_fits = fits.open(files[args.index])[0]
#event_mjd = skymap_fits.header['TIME_MJD']
event_mjd = 58000.000 #HARDCODE SO THAT THERE IS REAL DATA
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
    f = FastResponseAnalysis(files[args.index], start, stop, save=False, alert_event=True)
    inj = f.initialize_injector(gamma=gamma)
    for nsig in nsigs:
        for jj in range(trials_per_sig):
            seed_counter += 1
            ni, sample, true_ra, true_dec = inj.sample(nsig, poisson = True, return_position = True)
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
                flux_list.append(inj.mu2flux(nsig))
                true_ras.append(true_ra[0])
                true_decs.append(true_dec[0])

results = {'ts_prior': tsList_prior, 'ts': tsList, 'ns_prior': nsList_prior,
            'ns': nsList, 'true_ns': true_ns, 'ra': ra, 'dec': dec,
            'gamma': gammaList, 'mean_ninj': mean_ninj, 'flux': flux_list,
            'true_ra': true_ras, 'true_dec': true_decs}

with open(output_paths + 'index_{}_time_{}.pkl'.format(args.index, args.deltaT), 'w') as fi:
    pickle.dump(results, fi, protocol=pickle.HIGHEST_PROTOCOL)

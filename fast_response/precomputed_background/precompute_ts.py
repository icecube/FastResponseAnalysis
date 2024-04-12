#!/usr/bin/env python

r"""
Run trials for background only, all-sky scans. Record TS and 
number of true and fitted events around best fit location

"""
import numpy as np
import healpy as hp
import os, sys, argparse
from astropy.time import Time
from numpy.lib.recfunctions import append_fields
from scipy import sparse
from glob import glob

from fast_response.AlertFollowup import AlertFollowup

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--deltaT', type=float, default=None,
                    help='Time Window in seconds')
parser.add_argument('--ntrials', type=int, default = 10000,
                        help='Trials')
parser.add_argument("--bkg", default=6.4, type=float,
                        help="Expected background rate in mHz (default 6.4)")
parser.add_argument('--seed', default=1, type=int, 
                        help='Unique seed for running on the cluster')
parser.add_argument('--outdir',type=str, default=None,
                        help='Output directory to save npz (default = FAST_RESPONSE_OUTPUT env variable or cwd)')
args = parser.parse_args()

if args.outdir == None:
    try: 
        outdir = os.environ.get('FAST_RESPONSE_OUTPUT')
    except:
        outdir = './'
else:
    outdir = args.outdir
if not os.path.exists(outdir+'/trials/'):
    os.mkdir(outdir+'/trials/')
outdir=outdir+'/trials/'

skymap_files = glob('/data/ana/realtime/alert_catalog_v2/2yr_prelim/fits_files/Run13*.fits.gz')

#start_mjd = 58484.0 #Jan 1, 2019
start_mjd = 60298.0 #Dec 20, 2023
stop_mjd = start_mjd + (args.deltaT / 86400.)
start_iso = Time(start_mjd, format='mjd').iso
stop_iso = Time(stop_mjd, format='mjd').iso
deltaT = args.deltaT / 86400.

#skymap required for initialization, but not used here
f = AlertFollowup('Precompute_trials_test', skymap_files[0],
                start_iso, stop_iso, save=False)  
f.llh.nbackground=args.bkg*args.deltaT/1000.
#inj = f.initialize_injector(gamma=2.5) #just put this here to initialize f.spatial_prior
#print f.llh.nbackground
#results_array = []

ntrials = args.ntrials
stop = ntrials * (args.seed+1)
start = stop-ntrials
seed_counter = start

npix = hp.nside2npix(f.nside)
shape = (args.ntrials, npix)
maps = sparse.lil_matrix(shape, dtype=float)
for jj in range(ntrials):
    seed_counter += 1
    val = f.llh.scan(0.0, 0.0, scramble=True, seed = seed_counter,
            #spatial_prior = f.spatial_prior, 
            time_mask = [deltaT / 2., (start_mjd + stop_mjd) / 2.],
            pixel_scan = [f.nside, f._pixel_scan_nsigma], inject = None)
    if val['TS'] is not None:
        dtype = [('ts',float),('pixel',float)]
        results = np.empty((val['TS'].size,), dtype=dtype)
        pixels = hp.ang2pix(f.nside, np.pi/2. - val['dec'], val['ra'])
        maps[jj, pixels] = val['TS']
print("DONE")
hp_sparse = maps.tocsr()

outfilename = outdir+'{:.1f}_mHz_seed_{}_delta_t_{:.1e}.npz'.format(args.bkg, args.seed, args.deltaT)
sparse.save_npz(outfilename, hp_sparse)
print("Saved to {}".format(outfilename))
r"""
Run trials for background only, all-sky scans. Record TS and 
number of true and fitted events around best fit location

"""
import numpy as np
import healpy as hp
import os, sys, argparse
from astropy.time import Time
from astropy.time import TimeDelta
from numpy.lib.recfunctions import append_fields
from astropy.io import fits
from glob import glob
from scipy import sparse
import pickle

base_path = os.path.join('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/','')
sys.path.append(base_path)
from FastResponseAnalysis import FastResponseAnalysis

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--deltaT', type=float, default=None,
                    help='Time Window in seconds')
parser.add_argument('--ntrials', type=int, default = 10000,
                        help='Trials')
parser.add_argument("--bkg", default=6.4, type=float,
                        help="Expected background rate in mHz (defualt 6.4)")
parser.add_argument('--seed', default=1, type=int, 
                        help='Unique seed for running on the cluster')
args = parser.parse_args()

skymap_files = glob('/data/ana/realtime/alert_catalog_v2/2yr_prelim/fits_files/Run13*.fits.gz')

start_mjd = 58484.0
stop_mjd = start_mjd + (args.deltaT / 86400.)
start = Time(start_mjd, format='mjd').iso
stop = Time(stop_mjd, format='mjd').iso
deltaT = args.deltaT / 86400.

trials_per_sig = args.ntrials
seed_counter = args.seed

f = FastResponseAnalysis(skymap_files[0], start, stop, save=False,
                            alert_event=True)
f.llh.nbackground=args.bkg*args.deltaT/1000.
#inj = f.initialize_injector(gamma=2.5) #just put this here to initialize f.spatial_prior
#print f.llh.nbackground

#results_array = []

npix = hp.nside2npix(f.nside)
shape = (args.ntrials, npix)
maps = sparse.lil_matrix(shape, dtype=float)
import time
t0 = time.time()
for jj in range(trials_per_sig):
    seed_counter += 1
    val = f.llh.scan(0.0, 0.0, scramble=True, seed = seed_counter,
            #spatial_prior = f.spatial_prior, 
            time_mask = [deltaT / 2., (start_mjd + stop_mjd) / 2.],
            pixel_scan = [f.nside, 4.0], inject = None)
    if val['TS'] is not None:
        dtype = [('ts',float),('pixel',float)]
        results = np.empty((val['TS'].size,), dtype=dtype)
        pixels = hp.ang2pix(f.nside, np.pi/2. - val['dec'], val['ra'])
        maps[jj, pixels] = val['TS']
t1 = time.time()
print("took {:1f} seconds to do {} trials".format(t1-t0, trials_per_sig))
print("DONE")
hp_sparse = maps.tocsr()
outfilename = '/data/user/apizzuto/fast_response_skylab/fast-response/trunk/precomputed_background/trials/{:.1f}_mHz_seed_{}_delta_t_{:.1e}.npz'.format(args.bkg, args.seed, args.deltaT)
sparse.save_npz(outfilename, hp_sparse)
print("Saved to {}".format(outfilename))

#np.save('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/precomputed_background/trials/{:.1f}_mHz_seed_{}_delta_t_{:.1e}.npy'.format(args.bkg, args.seed, args.deltaT),
#        results_array)

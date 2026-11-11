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

# PS LLH without prior works for these BG scans
from fast_response.MultiExternalFollowup import MultiFollowup, GFUFollowup, GrecoFollowup, DNNFollowup, DNNOnlineFollowup, DNNIceManFollowup

parser = argparse.ArgumentParser(description='Precompute MultiAlertFollowup BG trials')
parser.add_argument('--deltaT', type=float, default=1000.,
                    help='Time Window in seconds')
parser.add_argument('--ntrials', type=int, default = 1000,
                        help='Trials')
parser.add_argument('--start', type=str, required=False,
                    default='2025-06-01',
                    help="Start time of the analysis in ISO format")
parser.add_argument("--bkg", default=[0.206], type=float, nargs='+',
                        help="Expected background rates in mHz (default 6.4, 4.6)")
parser.add_argument('--seed', default=1, type=int, 
                        help='Unique seed for running on the cluster')
parser.add_argument('--nside', default=128, type=int, 
                        help='Skymap nside to scan')
parser.add_argument('--outdir',type=str, default=os.environ.get('FAST_RESPONSE_OUTPUT', './'),
                        help='Output directory to save npz (default = FAST_RESPONSE_OUTPUT env variable or cwd)')
parser.add_argument('--dataset', default="DNN", type=str,
                        help='Dataset(s) to include, joined by + signs')
parser.add_argument('--fix_index', action='store_true', help='Fix the spectral index during fitting')
parser.add_argument('--index', type=float, default=None,
                        help="Spectral index to assume in injected hypothesis")
args = parser.parse_args()

# set up combination of datasets
followups = []
if 'GFU' in args.dataset:
    followups.append(GFUFollowup)
if 'Greco' in args.dataset:
    followups.append(GrecoFollowup)
if 'DNN' in args.dataset:
    followups.append(DNNOnlineFollowup)
if 'IceMan' in args.dataset:
    followups.append(DNNIceManFollowup)
MultiFollowup._followups = followups

# set fixed/floating index
for attr in ['index', 'fix_index']:
    setattr(MultiFollowup, f'_{attr}', getattr(args, attr))
MultiFollowup._float_index = not MultiFollowup._fix_index


# use the signal window to select events, as in a real analysis
start_mjd = Time(args.start).mjd
stop_mjd = start_mjd + (args.deltaT / 86400.)
start_iso = args.start
stop_iso = Time(stop_mjd, format='mjd').iso
deltaT = args.deltaT / 86400.

#skymap required for initialization, but not used here
f = MultiFollowup('Precompute_trials_test', 0,0,
                start_iso, stop_iso, save=False)
assert f.scramble
f._ncpu = 2

outdir = os.path.join(args.outdir, "precomputed_trials/", "+".join(f.datasets))
if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

# set the BG expectation assumed by the LLH
# NOTE the BG rate represented by the scrambles comes from args.start
for i, enum in enumerate(f.llh._samples):
    _llh = f.llh._samples[enum]
    # set background to required rate [mHz] to simulate
    _llh.nbackground = (args.bkg[i]/1000.) * _llh.on_livetime*86400
    print(_llh.nbackground)

# per original seed of the job, carve out a block of unique seeds, one per trial
# this permits reproducing these trials later
seeds = range(args.seed, args.seed + args.ntrials)

# set up skymaps
npix = hp.nside2npix(args.nside)
shape = (args.ntrials, npix)
maps = sparse.lil_matrix(shape, dtype=float)
for jj, seed in enumerate(seeds):
    val = f.llh.scan(0.0, 0.0, scramble=True, seed = seed,
            time_mask = [deltaT / 2., (start_mjd + stop_mjd) / 2.],
            pixel_scan = [args.nside, 3.0], inject = None)
    if val['TS'] is not None: 
        dtype = [('ts',float),('pixel',float)]
        results = np.empty((val['TS'].size,), dtype=dtype)
        pixels = hp.ang2pix(args.nside, np.pi/2. - val['dec'], val['ra'])
        maps[jj, pixels] = val['TS']
    else: # TODO same as in other scripts, but why would this happen?
        raise ValueError("scan TS is None")
print("DONE")
hp_sparse = maps.tocsr()

rates_str = '_'.join([f'{_rate:.2f}' for _rate in args.bkg])
outfilename = '_'.join([
    f'precomputed_trials_delta_t_{args.deltaT:.2e}',
    f'nside_{args.nside}',
    f'index_{f._index}',
    f'{rates_str}_mHz',
    f'seed_{args.seed}',
    ]
)
outfilepath = os.path.join(outdir, outfilename)
sparse.save_npz(outfilepath, hp_sparse)
print("Saved to {}".format(outfilepath))

#!/usr/bin/env python

r"""
Run trials for background only, all-sky scans. Record TS and 
number of true and fitted events around best fit location.
Datasets and their combinations can be specified.

Author: Christoph Raab
Date: 2026
"""
import numpy as np
import healpy as hp
import os, sys, argparse
from astropy.time import Time, TimeDelta
from numpy.lib.recfunctions import append_fields
from scipy import sparse
from glob import glob

# PS LLH without prior works for these BG scans
from fast_response.MultiExternalFollowup import MultiFollowup, GFUFollowup, GrecoFollowup, DNNFollowup, DNNOnlineFollowup, DNNIceManFollowup

parser = argparse.ArgumentParser(description='Precompute MultiAlertFollowup BG trials')
parser.add_argument('--deltaT', type=float, default=1000.,
                    help='Time Window [seconds]')
parser.add_argument('--ntrials', type=int, default = 1000,
                        help='Trials')
parser.add_argument('--time', type=str, required=False, default="2026-06-30T00:00:00",
                    help="GW time [ISO-T]")
parser.add_argument("--bkg", default=None, type=float, nargs='+',
                        help="Space-separated list of expected background rate per dataset in mHz (default: use sample estimates")
parser.add_argument('--seed', default=np.random.randint(2**32), type=int, 
                        help='Unique seed for running on the cluster')
parser.add_argument('--nside', default=256, type=int, 
                        help='Skymap nside to scan')
parser.add_argument('--outdir',type=str, default=os.environ.get('FAST_RESPONSE_OUTPUT', './'),
                        help='Output directory to save npz (default = FAST_RESPONSE_OUTPUT env variable or cwd)')
parser.add_argument('--dataset', default="DNN", type=str,
                        help='Dataset(s) to include, joined by + signs')
parser.add_argument('--fix_index', action='store_true', help='Fix the spectral index during fitting')
parser.add_argument('--index', type=float, default=2.0,
                        help="Spectral index to assume in injected hypothesis")

def main(args):

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
    time = Time(args.time, format="isot")
    delta_t = TimeDelta(args.deltaT, format="sec")
    start_time = time - delta_t / 2
    stop_time = time + delta_t / 2
    start_iso, stop_iso = start_time.iso, stop_time.iso
    start_mjd, stop_mjd = start_time.mjd, stop_time.mjd
    

    #skymap required for initialization, but not used here
    f = MultiFollowup('Precompute_trials_test', 0,0,
                    start_iso, stop_iso, save=False)
    assert f.scramble
    f._ncpu = 2

    outdir = os.path.join(args.outdir, "precomputed_trials/", "+".join(f.datasets))
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    # set the BG expectation assumed by the LLH
    # NOTE the BG rate represented by the scrambles comes from args.time
    if args.bkg is not None:
        # passed an argument
        for i, enum in enumerate(f.llh._samples):
            _llh = f.llh._samples[enum]
            # set background to required rate [mHz]
            _llh.nbackground = (args.bkg[i] / 1000.) * (_llh.on_livetime * 86400)
    else:
        # other direction, set from samples [mHz]
        # which set nbackground = bkg_rate * on_livetime
        args.bkg = [1000 * _llh.nbackground / (_llh.on_livetime * 86400)
                    for enum, _llh in f.llh._samples.items()]
        # => will be used to label files

    # per original seed of the job, carve out a block of unique seeds, one per trial
    # this permits reproducing these trials later
    seeds = range(args.seed, args.seed + args.ntrials)

    # set up skymaps
    npix = hp.nside2npix(args.nside)
    shape = (args.ntrials, npix)
    # since many pixels don't get a value, sparse matrix is most space efficient
    maps = sparse.lil_matrix(shape, dtype=float)
    for jj, seed in enumerate(seeds):
        val = f.llh.scan(0.0, 0.0, scramble=True, seed = seed,
                time_mask = [(stop_mjd - start_mjd) / 2., (start_mjd + stop_mjd) / 2.],
                pixel_scan = [args.nside, 3.0], inject = None)
        if val['TS'] is not None: 
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
    print(f"Saved to {outfilepath}")

if __name__ == "__main__":
    args = parser.parse_args()
    sys.exit(main(args))
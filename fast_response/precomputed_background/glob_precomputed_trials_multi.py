#!/usr/bin/env python
"""
This module allows combining precomputed background sky scans
from multiple files saves as SciPy sparse matrices.
This is based on glob_precomputed_trials.py but generalized
in the assumptions on the file name pattern, so it can be used 
with different conventions, the only requirement that they contain
some seed index in the form `seed_*`.
Running as a script, it will save the combined maps again as a sparse matrix.
"""

import logging
from glob import glob
import healpy as hp
from scipy import sparse
import time
import argparse
import os
import sys

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description='Glob precomputed trials')
parser.add_argument('--dir',type=str, default='./',
                    help='directory for where trials are, will save globbed npz inside same')
parser.add_argument('--nside',type=int, default=256,
                    help='nside used when running trials (default 256)')

def get_glob_file(filename):
    dirname = os.path.dirname(filename)
    basename = os.path.basename(filename)
    particles = basename.replace('seed_','seed').split('_')
    basename_glob = '_'.join(_p for _p in particles if not 'seed' in _p)
    return os.path.join(dirname, basename_glob)

def sort_by_glob_file(files):
    sorted = {}
    for _file in files:
        outfile = get_glob_file(_file)
        sorted.setdefault(outfile, [])
        sorted[outfile].append(_file)
    return sorted

def load_maps(fn):
    return sparse.load_npz(fn)

def concatenate_maps(files, nside) -> sparse.csr_matrix:
    files = [fn for fn in files if f"nside_{nside}" in fn]
    logger.info('Found {} files to load'.format(len(files)))
    if len(files)==0:
        return None
    logger.info('Nside: {}'.format(nside))
    npix = hp.nside2npix(nside)
    logger.info('Starting to load at {}'.format(time.ctime()))    
    maps = sparse.csr_matrix((0, npix), dtype=float)
    for fn in files:
        scan = load_maps(fn)
        maps = sparse.vstack((maps, scan))
    return maps

def save_maps(maps, out):
    logger.info(f'Creating {os.path.basename(out)}')
    # Change format of sparse array
    logger.info("Starting to change from COO to CSR at {}".format(time.ctime()))
    # NOTE: I am not sure we need to keep this; maps are already CSR
    scans = maps.tocsr()
    logger.info("Finished at {}".format(time.ctime()))
    # Save the sparse array
    sparse.save_npz(out, scans)


def main(args):
    """
    Glob the all-sky scans together
    """
    # separate alert and GW trials by directory
    files = sorted(glob(os.path.join(args.dir, '*seed_*.npz')))
    for outfile, filegroup in sort_by_glob_file(files).items():
        maps = concatenate_maps(filegroup, args.nside)
        save_maps(maps, outfile)
        del maps
    logger.info('done')
    return 0

if __name__ == "__main__":
    args = parser.parse_args()
    logger.setLevel(logging.INFO)
    sys.exit(main(args))
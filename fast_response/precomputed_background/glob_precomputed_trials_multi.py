#!/usr/bin/env python

from glob import glob
import healpy as hp
from scipy import sparse
import time
import argparse
import os

parser = argparse.ArgumentParser(description='Glob precomp trials')
parser.add_argument('--dir',type=str, default='./',
                    help='directory for where trials are, will save globbed npz inside same ')
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


    
def glob_allsky_scans(files, out, nside):
    print(f'Creating {os.path.basename(out)}')
    files = [fn for fn in files if f"nside_{nside}" in fn]
    print('Found {} files to load'.format(len(files)))
    if len(files)==0: return None
    print('Nside: {}'.format(nside))
    npix = hp.nside2npix(nside)
    print('Starting to load at {}'.format(time.ctime()))    
    maps = sparse.csr_matrix((0, npix), dtype=float)
    for f in files:
        scan = sparse.load_npz(f)
        maps = sparse.vstack((maps, scan))
        print('.', end=' ')
    print('')

    # Change format of sparse array
    print("Starting to change from COO to CSR at {}".format(time.ctime()))
    scans = maps.tocsr()
    print("Finished at {}".format(time.ctime()))
    
    # Save the sparse array
    sparse.save_npz(out, scans)
    return maps

def main(args):
    """
    Glob the all-sky scans together
    """
    # separate alert and GW trials by directory
    files = sorted(glob(os.path.join(args.dir, '*seed_*.npz')))
    for outfile, filegroup in sort_by_glob_file(files).items():
        maps = glob_allsky_scans(filegroup, outfile, args.nside)
        del maps
    print ('done')

args = parser.parse_args()
main(args)
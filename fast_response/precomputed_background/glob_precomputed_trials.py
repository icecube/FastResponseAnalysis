#!/usr/bin/env python

from glob import glob
import healpy as hp
from scipy import sparse
import numpy as np
import time
import pickle
import argparse
import os

parser = argparse.ArgumentParser(description='Glob precomp trials')
parser.add_argument('--deltaT', type=float, default=None,
                    help='Time Window in seconds')
parser.add_argument('--dir',type=str, default='./',
                    help='directory for where trials are, will save npz to dir+/glob_trials/')
parser.add_argument('--nside',type=int, default=256,
                    help='nside used when running trials (default 256)')
parser.add_argument('--type', type=str, default='alert',
                    help='type of trials run (alert or gw)')
args = parser.parse_args()

def glob_allsky_scans(delta_t, rate, dir, low_stats=False, typ='alert'):
    """
    Glob the all-sky scans together
    """
    #n_per_job = {1000.: 1000, 172800.: 50, 2678400.: 30}
    #jobs_per_window = {1000.: 20, 172800.: 100, 2678400.: 100}

    if typ == 'gw':
        files = glob(dir+'/gw_{:.1f}_mHz_delta_t_{:.1e}_seed_*.npz'.format(rate, delta_t))
    elif typ =='alert':
        files = glob(dir+'/{:.1f}_mHz_seed_*_delta_t_{:.1e}.npz'.format(rate, delta_t))
    else:
        print('Error: type is not alert or gw')
        return
    
    nside = args.nside
    print('Nside: {}'.format(args.nside))
    npix = hp.nside2npix(nside)
    maps = sparse.csr_matrix((0, npix), dtype=float)
    
    print('Found {} files to load'.format(len(files)))
    if len(files)==0: return None
    
    print('Starting to load at {}'.format(time.ctime()))
    final_ind = -1 if not low_stats else 10
    for f in files[:final_ind]:
        scan = sparse.load_npz(f)
        maps = sparse.vstack((maps, scan))

    # Change format of sparse array
    print("Starting to change from COO to CSR at {}".format(time.ctime()))
    scans = maps.tocsr()
    print("Finished at {}".format(time.ctime()))
    
    # Save the sparse array
    stats_str = '' if not low_stats else '_low_stats'
    outfilename = '{}_precomputed_trials_delta_t_{:.2e}_trials_rate_{:.1f}{}.npz'.format(
        typ, delta_t, rate, stats_str)
    save_dir = os.path.join(dir, 'glob_trials')
    
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    #scan_path = '/data/user/apizzuto/fast_response_skylab/fast-response/fast_response/precomputed_background/glob_trials/'
    sparse.save_npz(os.path.join(save_dir,outfilename), scans)
    
    return maps

for rate in [6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2]:
    for low_stats in [True, False]:
        print("Rate: {} mHz, low stats: {}".format(rate, low_stats))
        maps = glob_allsky_scans(args.deltaT, rate, args.dir, typ=args.type, low_stats=low_stats)
        del maps
        print ('done')


#for rate in [6.2, 6.4, 6.6, 6.8, 7.0, 7.2]:
#    for low_stats in [True, False]:
#        print("Rate: {} mHz, low stats: {}".format(rate, low_stats))
#        maps = glob_allsky_scans(1000., rate, low_stats=low_stats)
#        del maps
#        print ('')

#for rate in [7.0, 7.2]:
#    for low_stats in [True]:
#        print("Rate: {} mHz, low stats: {}".format(rate, low_stats))
#        maps = glob_allsky_scans(2.*86400., rate, low_stats=low_stats)
#        del maps
#        print ('')

#!/usr/bin/env python

'''
Script to generate ntrials number of background only test statistics and save 
information about each trial to outfile

@options: 
    --i index: alert event index
    --ntrials: number of trials to perform
'''
import os, time, sys, pickle, argparse
import numpy as np
import healpy as hp

from skylab.ps_llh              import PointSourceLLH, MultiPointSourceLLH
from skylab.ps_injector         import PriorInjector
from skylab.llh_models          import EnergyLLH
from skylab.datasets            import Datasets
from skylab.priors              import SpatialPrior
from astropy.io                 import fits

sys.path.append('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/time_integrated_scripts/')
from config_steady              import config

##################################### CONFIGURE ARGUMENTS #############################
parser = argparse.ArgumentParser(description = 'Alert event followup steady background')
parser.add_argument('--i', type=int, required=True, help='Alert event index')
parser.add_argument('--verbose', action='store_true', default=False,
                    help="Assorted print statements flag")
parser.add_argument('--rng', type=int, default=1, help="Random number seed")
parser.add_argument('--smear', default=False, action='store_true',
                    help='Include systematics by smearing norm. prob.')
parser.add_argument('--local_skymap', default=False, action='store_true',
                    help='Instead of returning best fit, return whole map')
args = parser.parse_args()
#######################################################################################

index = args.i
seed = args.rng
verbose = args.verbose

smear_str = 'smeared/' if args.smear else 'norm_prob/'
if not args.local_skymap:
    outfile = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/results/{}index_{}_steady_seed_{}.pkl'.format(smear_str, index, seed)
else:
    outfile = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/results/{}index_{}_steady_ts_map_seed_{}.pkl'.format(smear_str, index, seed)

t0 = time.time()

nside = 2**7
multillh, spatial_prior = config(index, gamma = 2.0, seed = seed, scramble = False, nside=nside, 
                        ncpu = 1, injector = False, verbose=verbose, smear=args.smear)

t1 = time.time()
if verbose:
    print('{:.2f} seconds to Initialize Likelihoods'.format(t1 - t0))
    print ("\nRunning fit on real data ...")

allspots = None
ii = 1
n_iter = 2 if not args.local_skymap else 1
for results, hotspots in multillh.do_allsky_trials(n_iter= n_iter, 
                              injector=None,
                              nside=nside, rng_seed = 123*seed + ii,
                              spatial_prior=spatial_prior,
                              follow_up_factor = 1,
                              scramble = False):
    if verbose:
        print('Trial Number: {}'.format(ii))
    ii += 1
    if not args.local_skymap:
        if allspots is None:
            allspots = {}
            for k, v in hotspots['spatial_prior_0']['best'].items():
                allspots[k] = [v]
            if 'pix' not in allspots.keys():
                allspots['pix'] = [0]
            if 'nside' not in allspots.keys():
                allspots['nside'] = [0]
        else:
            for k, v in hotspots['spatial_prior_0']['best'].items():
                allspots[k].append(v)
            if 'pix' not in hotspots['spatial_prior_0']['best'].keys():
                allspots['pix'].append(0)
            if 'nside' not in hotspots['spatial_prior_0']['best'].keys():
                allspots['nside'].append(0)
    else:
        allspots = results
    #allspots.append(hotspots)

if args.local_skymap:
    allspots = allspots[allspots['TS'] != 0.]

dt1 = t1 - t0
dt     = time.time() - t0
if verbose:
    print("Finished script in {} seconds".format(dt))
    print("Initialization: {} seconds\ntrials: {} seconds".format(dt1, (dt-dt1)))

with open(outfile, 'w') as f:
    pickle.dump(allspots, f, protocol=pickle.HIGHEST_PROTOCOL)

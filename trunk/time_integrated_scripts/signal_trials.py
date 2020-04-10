#!/usr/bin/env python

'''
Signal injection for time integrated ANITA analysis

@options:
    --i:    rng_seed
    --g:    Spectral Index
    --ntrials: number of trials
    --nsignal: mean injected signal in number of events
    --event: ANITA event
'''

import os, logging, inspect, argparse, time, sys, pickle

from scipy.stats import chi2
import healpy as hp
import numpy as np
from scipy.special import erfc
from glob import glob


from skylab.ps_llh              import PointSourceLLH, MultiPointSourceLLH
from skylab.ps_injector         import PointSourceInjector, PriorInjector
from skylab.llh_models          import EnergyLLH, ClassicLLH
from skylab.datasets            import Datasets
from skylab.priors              import SpatialPrior
from skylab.sensitivity_utils   import fit
from skylab.test_statistics     import TestStatisticNotSeparated

sys.path.append('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/time_integrated_scripts/')
from config_steady              import config


###############################################################################
parser = argparse.ArgumentParser(description = 'ANITA_Tau_Sens_seed_timing_tests')
#parser.add_argument('--nsignal', type=float, required=True, help='mean number of injected signal events')
parser.add_argument('--i', type=int, required=True, help='Alert event index')
parser.add_argument('--g', type=float, required=True, help='spectral index')
parser.add_argument("--ntrials", default=10, type=int,
                help="Number of trials per signal strength (default=10)")
parser.add_argument('--rng', type=int, default=1, help="Random number seed")
parser.add_argument('--verbose', action='store_true', default=False,
                    help="Assorted print statements flag")   
parser.add_argument('--fit', action='store_true', default=False,
                    help="Include poisson fluctuations by default or raise this flag")                
args = parser.parse_args()
###############################################################################

gamma = args.g
index = args.i
ntrials = args.ntrials
seed = args.rng
verbose = args.verbose
poisson = ~args.fit

outdir = 'fits' if args.fit else 'sensitivity'
outfile = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/{}/index_{}_steady_seed_{}_gamma_{}.pkl'.format(outdir, index, seed, gamma)

t0 = time.time()

nside = 2**7
multillh, spatial_prior, inj = config(index, gamma = gamma, seed = seed, scramble = True, nside=nside,
                        ncpu = 1, injector = True, verbose=verbose)


t1 = time.time()
if verbose:
    print('{:.2f} seconds to Initialize Likelihoods'.format(t1 - t0))
    print ("\nRunning signal injection trials ...")

allspots = None
ii = 1
ninj = np.append(np.linspace(1, 10, 10), np.array([15, 20, 25, 30, 40, 50, 60, 70]))

scale_arr = []
for i in range(1,21):
    scale_arr.append([])
    for j in range(5):
        scale_arr[-1].append(inj.sample(i, poisson=False)[0][0])
scale_arr = np.median(scale_arr, axis=1)
try:
    scale_factor = np.min(np.argwhere(scale_arr > 0)) + 1.
except:
    print("Scale factor thing for prior injector didn't work")
    scale_factor = 1.

for ni in ninj:
    for results, hotspots in multillh.do_allsky_trials(n_iter= ntrials,
                                injector=inj,
                                mean_signal=ni*scale_factor,
                                poisson=poisson, 
                                nside=nside, rng_seed = 123*seed + ii,
                                spatial_prior=spatial_prior,
                                follow_up_factor = 1, return_position=True):
        if verbose:
            print('Trial Number: {}'.format(ii))
        ii += 1
        if allspots is None:
            allspots = {}
            for k, v in hotspots['spatial_prior_0']['best'].items():
                allspots[k] = [v]
            if 'pix' not in allspots.keys():
                allspots['pix'] = [0]
            if 'nside' not in allspots.keys():
                allspots['nside'] = [0]
            if 'inj' in hotspots['spatial_prior_0'].keys():
                for k, v in hotspots['spatial_prior_0']['inj'].items():
                    allspots['inj_' + k] = [v]
            else:
                allspots['inj_nsignal'] = [0]
                allspots['inj_dec'], allspots['inj_ra'] = [0.], [0.]
            allspots['flux'] = [inj.mu2flux(ni*scale_factor)]
        else:
            for k, v in hotspots['spatial_prior_0']['best'].items():
                allspots[k].append(v)
            if 'pix' not in hotspots['spatial_prior_0']['best'].keys():
                allspots['pix'].append(0)
            if 'nside' not in hotspots['spatial_prior_0']['best'].keys():
                allspots['nside'].append(0)
            if 'inj' in hotspots['spatial_prior_0'].keys():
                for k, v in hotspots['spatial_prior_0']['inj'].items():
                    allspots['inj_' + k].append(v)
            else:
                allspots['inj_nsignal'].append(0)
                allspots['inj_dec'].append(0.0)
                allspots['inj_ra'].append(0.0)
            allspots['flux'].append(inj.mu2flux(ni*scale_factor))
        #allspots.append(hotspots)

dt1 = t1 - t0
dt     = time.time() - t0
if verbose:
    print("Finished script in {} seconds".format(dt))
    print("Initialization: {} seconds\ntrials: {} seconds".format(dt1, (dt-dt1)))

with open(outfile, 'w') as f:
    pickle.dump(allspots, f, protocol=pickle.HIGHEST_PROTOCOL)

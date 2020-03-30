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

import os
import logging
import inspect

# scipy
from scipy.stats import chi2
import healpy as hp
import numpy as np
from scipy.special import erfc

# local
#import utils
#import data

from glob import glob
from time import time

from matplotlib import mlab

import argparse
from skylab.ps_llh              import PointSourceLLH, MultiPointSourceLLH
from skylab.ps_injector         import PointSourceInjector, PriorInjector
from skylab.llh_models          import EnergyLLH, ClassicLLH
from skylab.datasets            import Datasets
from skylab.priors              import SpatialPrior
from skylab.sensitivity_utils   import fit
from skylab.test_statistics     import TestStatisticNotSeparated

from config_steady_ANITA        import config, Askaryan, Tau, Tau_2


###############################################################################
parser = argparse.ArgumentParser(description = 'ANITA_Tau_Sens_seed_timing_tests')
parser.add_argument('--nsignal', type=float, required=True, help='mean number of injected signal events')
parser.add_argument('--i', type=int, required=True, help='rng_seed')
parser.add_argument('--g', type=float, required=True, help='spectral index')
parser.add_argument("--ntrials", default=1000, type=int,
                help="Number of trials (default=1000")
parser.add_argument("--event", default="Askaryan", type=str,
                help="ANITA event location (Askaryan or Tau)")
args = parser.parse_args()
###############################################################################

if args.event == "Askaryan":
    source = Askaryan()
elif args.event == "Tau":
    source = Tau()
elif args.event == "Tau_2":
    source = Tau_2()
else:
    print("Invalid ANITA Event.\nExiting")
    exit(0)

rng_seed = args.i
gamma = args.g
ns = args.nsignal
ntrials = args.ntrials

outfile = '/data/user/apizzuto/ANITA/SteadyTrials/signal_injection/{}/gamma_{}_ns_{}_sig_trial_{:03d}'.format(args.event, gamma, ns, rng_seed)


t0 = time()

GeV = 1
TeV = 1000 * GeV

#################
# SPATIAL PRIOR #
#################

nside = 2**6
npix = hp.nside2npix(nside)
pixels = range(npix)

spatial_prior = SpatialPrior(source['prior_map'], scan_thresh = source['threshold'])


##########
# SKYLAB #
##########

seasons = [ ["PointSourceTracks_v003p01", "IC86, 2011"],
            ["PointSourceTracks_v003p01", "IC86, 2012-2017"] ]

multillh, inj = config(seasons, gamma = gamma, seed = rng_seed, scramble = True,
                        ncpu = 1, spatial_prior = spatial_prior, injector = True)

print("\n injected spectrum:")
print("   - %s" % str(inj.spectrum))

t1 = time()
print('\n{:.2f} seconds to Initialize Likelihoods'.format(t1 - t0))

##########
# SKYLAB #
##########

###############################################################################

#################
# SIGNAL TRIALS #
#################

print("\nSignal Trials:")


# NOTE: mean_signal *is not* the number of injected events in the case
#       of the prior injector because the injector adjusts the number
#       of events depending on the declination of injection to account
#       for the detector's effective area.
#
#       prior injector should be thought of as a constant flux injector.
#       This is why we convert mean_signal to flux and save that instead
#       of the number of injected events.

allspots = []
for results, hotspots in multillh.do_allsky_trials(n_iter=ntrials,  # run xx randomized trials
                                              rng_seed=rng_seed,    # pass different seed for each set of trials
                                              injector=inj,         # injector
                                              mean_signal=ns,       # signal to inject
                                              return_position=True, # track injected event info in hotspots dict
                                              nside=nside,          # nside of map used for initial, coarse search
                                              follow_up_factor=2,   # each successive map uses (2^1) * nside
                                              spatial_prior=spatial_prior): # apply spatial prior when looking for hotspot
    allspots.append(hotspots)

dt1 = t1 - t0
dt     = time() - t0
print("Finished script in {} seconds".format(dt))
print("Initialization: {} seconds\ntrials: {} seconds".format(dt1, (dt-dt1)))

np.save(outfile, allspots)


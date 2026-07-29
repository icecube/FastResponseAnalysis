'''Script to calculate sensitivities of analyses defined
in the fast_response framework.

    Author: Christoph Raab
    Date: 2026
    '''

# from fast_response.MultiGWFollowup import IceManFollowup as MultiFollowup
from fast_response.MultiGWFollowup import DNNOnlineFollowup as MultiFollowup
import argparse
import subprocess
import warnings
import fast_response.web_utils as web_utils
import logging as log
import pyfiglet
import os
import sys
import numpy as np
import scipy
import pickle
from tqdm import tqdm

from skylab import utils
from skylab.llh_models      import EnergyLLH
from skylab.ps_llh          import PointSourceLLH
from skylab.temporal_models import BoxProfile, TemporalModel

from astropy.time import Time, TimeDelta
import datetime


def run_trials(args):
    result_labels = ['prior_sensitivity_scan']
    outdir = {}
    skymap_name = os.path.split(args.skymap)[1]
    for _r in result_labels:
        _outdir = os.path.join(args.out,_r, skymap_name)
        if not os.path.exists(_outdir):
            os.makedirs(_outdir)
        outdir[_r] = _outdir
    

    for attr in ['index', 'fix_index']:
        setattr(MultiFollowup, f'_{attr}', getattr(args, attr))
    MultiFollowup._float_index = not MultiFollowup._fix_index
    
    llh_seed = np.random.randint(0, 1000000)
    f = MultiFollowup(
        args.name, args.skymap, args.start,
        args.stop, 
        skipped=args.skip_events, save=False,
    )
    assert f.scramble

    dataset_string = '+'.join(sorted(f.datasets))
    outpath = os.path.join(outdir["prior_sensitivity_scan"],
                           f"mu_{args.mu:.2f}_{dataset_string}.parquet")
    if os.path.exists(outpath):
        sys.exit()

    f._allow_neg = False
    # f._ncpu = args.ncpu
    # llh=config_llh(f)
    nside=f.nside

    print('Init with values:')
    print(f'_allow_neg == {f._allow_neg}')
    print(f'_ncpu == {f._ncpu}')
    print(f'nside == {f.nside}')

    print(str(f))
    print(f.dataset)
    print(args.out)
    print(f.outdir)
    print(f.save_output,
          f.analysisid,
          f.analysispath,
         )
    
    f.initialize_injector()
    spatial_prior = f.spatial_prior
    inj=f.inj
    s= np.random.randint(0, 100000)
    # different seed?
    inj.set_rng_seed(s)
    print('Initialized injector')

    ns = args.mu
    flux = inj.mu2flux(ns)

    TS_list, gamma_fit=[], []
    ns_fit, ns_inj =[], []
    flux_inj, flux_fit =[], []
    for j in tqdm(range(args.n_iter)):
        ni, sample = inj.sample(ns,poisson=True)
        val = f.llh.scan(0.0,0.0, scramble = True, seed = j, spatial_prior=spatial_prior,
                    inject = sample,time_mask=[f.duration/2., f.centertime], pixel_scan=[nside,3.])

        if val.size > 0:
            maxLoc = np.argmax(val['TS_spatial_prior_0'])  #pick out max of all likelihood ratios at diff pixels
            TS_list.append(val['TS_spatial_prior_0'].max())
        else:
            # TS_list.append(-np.inf)
            TS_list.append(0) # new convention
        
        
        if sample is not None: 
            sample = list(sample.values())[0]
            ns_inj.append(sample.size)
            flux_inj.append(inj.mu2flux(sample.size))
        else: 
            ns_inj.append(0.)
            flux_inj.append(0.)
        
        if val.size > 0:
            ns_fit.append(val['nsignal'][maxLoc])
            if f._float_index:
                gamma_fit.append(val['gamma'][maxLoc])
            flux_fit.append(inj.mu2flux(val['nsignal'][maxLoc]))
        else:
            ns_fit.append(0)
            flux_fit.append(0)
            if f._float_index:
                gamma_fit.append(f.index)
    
    import pandas as pd
    results = pd.DataFrame(
        {
        'TS_List':TS_list,
        'ns_fit':ns_fit,
        'ns_inj':ns_inj,
        'gamma_fit':gamma_fit,
        'flux_inj':flux_inj,
        'flux_fit':flux_fit,
        })
    for k, v in {
        'mean_ns':ns,
        'mean_flux':flux,
        'llh_seed':llh_seed,
        'inj_seed':s,
        }.items():
        results[k] = v

    
    results.to_parquet(outpath)
    
    

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    log.basicConfig(level=log.ERROR)

    parser = argparse.ArgumentParser(description='Fast Response Analysis sensitivity')
    parser.add_argument('--name', type=str, default="test inherit sensitivity",
                        help='Name of the source (do not use underscores or LaTeX will be mad!)')
    parser.add_argument("--skymap")
    parser.add_argument("--alpha", default=0.5, type=float)
    parser.add_argument("--beta", default=0.9, type=float)
    parser.add_argument('--fix_index', action='store_true', help='Fix the spectral index during fitting')
    parser.add_argument('--duration', default=10, type=float,
                        help='Duration of followup [days]')
    parser.add_argument('--start', type=str, required=False,
                        default=Time(58933.0 - 100, format='mjd').iso,
                        help="Start time of the analysis in ISO format")
    parser.add_argument('--extension', type=float, default=None,
                        help="Source extension in degrees")
    parser.add_argument('--index', type=float, default=None,
                        help="Spectral index to assume in injected hypothesis")
    
    parser.add_argument('--skip-events', default=None,
                        type= lambda z:[ tuple(int(y) for y in x.split(':')) for x in z.split(',')],
                        help="Event to exclude from the analyses, eg."
                        "Example --skip-events=127853:67093193")
    parser.add_argument('--ntrials', default=100000, type=int,
                        help="Number of background trials to perform")
    parser.add_argument('--n_iter', default=1000, type=int,
                        help="Number of signal trials per per sensitivity iteration")
    parser.add_argument('--out', default='/data/user/chraab/Output/fast_response/multisample/inherit/')
    parser.add_argument("--mu", type=float)
    
    args = parser.parse_args()

    if '_' in args.name:
        print('Warning: underscores in source name cause LaTeX to fail!')
    tic = datetime.datetime.now()

    args.stop = (Time(args.start, format='iso') + TimeDelta(args.duration, format='jd')).iso
    
    run_trials(args)
    toc = datetime.datetime.now()
    print((toc - tic).total_seconds())

'''Script to calculate sensitivities of analyses defined
in the fast_response framework.

    Author: Christoph Raab
    Date: 2026
    '''

# from fast_response.MultiGWFollowup import IceManFollowup as MultiFollowup
from fast_response.MultiGWFollowup import DNNOnlineFollowup as MultiFollowup
from fast_response.sensitivity_utils import binomial_error
import argparse
import subprocess
import warnings
import fast_response.web_utils as web_utils
import logging as log
import pyfiglet
import os
import pandas as pd
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
    result_labels = ['unblinding_scrambled']
    outdir = {}
    skymap_name = os.path.split(args.skymap)[1]
    for _r in result_labels:
        if args.logE is None:
            _outdir = os.path.join(args.outdir,_r, skymap_name)
        else:
            _outdir = os.path.join(args.outdir, _r, f"logE_{args.logE:.0f}", skymap_name)
        if not os.path.exists(_outdir):
            os.makedirs(_outdir)
        outdir[_r] = _outdir
    

    for attr in ['index', 'fix_index']:
        setattr(MultiFollowup, f'_{attr}', getattr(args, attr))
    MultiFollowup._float_index = not MultiFollowup._fix_index
    
    f = MultiFollowup(# from fast_response.MultiGWFollowup import IceManFollowup as MultiFollowup
        args.name, args.skymap, args.start,
        args.stop, 
        skipped=args.skip_events, save=False,
        seed=args.seed,
    )
    assert f.scramble

    f._allow_neg = False
    

    print('Init with values:')
    print(f'_allow_neg == {f._allow_neg}')
    print(f'_ncpu == {f._ncpu}')
    print(f'nside == {f.nside}')

    print(str(f))
    print(f.dataset)
    print(args.outdir)
    print(f.outdir)
    print(f.save_output,
          f.analysisid,
          f.analysispath,
         )
    

    # single-sample
    month = Time(args.start).datetime.month
    llh = f.llh._samples[0]
    # use the estimate provided by the temporal model
    # from the background window(s): (temporal model has estimated it from +/- 60 days)
    rate = round(1000 * llh.nbackground / (llh.on_livetime * 86400), 2)
    bg_ts = f.load_background_trials(rate=rate)

    results_list = []

    for seed in range(args.seed, args.seed + args.nseed):
        f.llh_seed = seed # used to scramble if scramble=True
        unblind_TS, unblind_ns, unblind_gamma = f.unblind_TS(scramble=True)

        results = {}
        
        results['rate'] = rate
        results['month'] = month
        results['duration'] = args.duration
        results['index'] = args.index
        results['ntrials'] = bg_ts.size
        results['seed'] = seed
        results['TS'] = unblind_TS,
        results['ns'] = unblind_ns,
        results['gamma'] = unblind_gamma,
        results['p'] = np.count_nonzero(bg_ts >= unblind_TS) / bg_ts.size,
        results['skymap'] = args.skymap,
        results_list.append(results)


    dataset_string = '+'.join(sorted(f.datasets))
    label = f"duration_{args.duration}_index_{args.index}_{dataset_string}_{rate:.2f}mHz_seed_{args.seed}"
    pd.DataFrame(results_list).to_parquet(
        os.path.join(outdir["unblinding_scrambled"],
                     label + ".parquet",
                           ))
    

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    log.basicConfig(level=log.ERROR)

    parser = argparse.ArgumentParser(description='Fast Response Analysis sensitivity')
    parser.add_argument('--name', type=str, default="test inherit sensitivity",
                        help='Name of the source (do not use underscores or LaTeX will be mad!)')
    parser.add_argument("--skymap", type=str, required=True)
    parser.add_argument('--fix_index', action='store_true', help='Fix the spectral index during fitting')
    parser.add_argument('--duration', default=1000./86400, type=float,
                        help='Duration of followup [days]')
    parser.add_argument('--start', type=str, required=False,
                        default='2026-06-30',
                        help="Start time of the analysis in ISO format")
    parser.add_argument('--index', type=float, default=2.0,
                        help="Spectral index to assume in injected hypothesis")
    parser.add_argument('--skip-events', default=None,
                        type= lambda z:[ tuple(int(y) for y in x.split(':')) for x in z.split(',')],
                        help="Event to exclude from the analyses, eg."
                        "Example --skip-events=127853:67093193")
    parser.add_argument('--ntrials', default=10000, type=int,
                        help="Number of background trials to load")
    parser.add_argument('--seed', default=100000, type=int,
                            help="RNG seed for the LLH")
    parser.add_argument("--nseed", type=int, default=1,
                        help="Number of seeds starting at --seed to use")
    parser.add_argument('--outdir', default='./')
    parser.add_argument("--logE", type=float, default=None,
                            help="set BG reco energy")
    
    args = parser.parse_args()

    if '_' in args.name:
        print('Warning: underscores in source name cause LaTeX to fail!')
    tic = datetime.datetime.now()

    args.stop = (Time(args.start, format='iso') + TimeDelta(args.duration, format='jd')).iso
    
    run_trials(args)
    toc = datetime.datetime.now()
    print((toc - tic).total_seconds())

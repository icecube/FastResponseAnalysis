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
    result_labels = ['prior_sensitivity_thresh', 'prior_bg_trials']
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
    
    f = MultiFollowup(# from fast_response.MultiGWFollowup import IceManFollowup as MultiFollowup
        args.name, args.skymap, args.start,
        args.stop, 
        skipped=args.skip_events, save=False,
    )
    assert f.scramble

    f._allow_neg = False
    

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

    month = Time(args.start).datetime.month
    # FIXME Horrible! do better than hardcoding this
    rate = [
        0.22, 0.22, 0.22, 0.21,
        0.21, 0.21, 0.21, 0.2 ,
        0.22, 0.23, 0.22, 0.23,
        ][month - 1]
    bg_ts = f.load_background_trials(rate=rate, ntrials=args.ntrials)

    results_list = []
    for alpha in args.alpha:
        results = {}
        results["TSval"] = np.sort(bg_ts)[int(bg_ts.size*(1-alpha))]
        results["alpha"] = alpha
        results["alpha_unc"] = binomial_error(alpha, bg_ts.size)
        results['rate'] = rate
        results['month'] = month
        results['duration'] = args.duration
        results['index'] = args.index
        results['ntrials'] = bg_ts.size
        results_list.append(results)


    dataset_string = '+'.join(sorted(f.datasets))
    label = f"duration_{args.duration}_index_{args.index}_{dataset_string}_{rate:.2f}mHz"
    pd.DataFrame(results_list).to_parquet(
        os.path.join(outdir["prior_sensitivity_thresh"],
                     label + ".parquet",
                           ))
    np.save(os.path.join(outdir["prior_bg_trials"], label + ".npy"), bg_ts)
    

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    log.basicConfig(level=log.ERROR)

    parser = argparse.ArgumentParser(description='Fast Response Analysis sensitivity')
    parser.add_argument('--name', type=str, default="test inherit sensitivity",
                        help='Name of the source (do not use underscores or LaTeX will be mad!)')
    parser.add_argument("--skymap")
    parser.add_argument("--alpha", default=[0.5, 0.00135], type=float, nargs='+')
    parser.add_argument('--fix_index', action='store_true', help='Fix the spectral index during fitting')
    parser.add_argument('--duration', default=1000./86400, type=float,
                        help='Duration of followup [days]')
    parser.add_argument('--start', type=str, required=False,
                        default='2025-06-01',
                        help="Start time of the analysis in ISO format")
    parser.add_argument('--index', type=float, default=2.0,
                        help="Spectral index to assume in injected hypothesis")
    parser.add_argument('--skip-events', default=None,
                        type= lambda z:[ tuple(int(y) for y in x.split(':')) for x in z.split(',')],
                        help="Event to exclude from the analyses, eg."
                        "Example --skip-events=127853:67093193")
    parser.add_argument('--ntrials', default=100000, type=int,
                        help="Number of background trials to load")
    parser.add_argument('--out', default='/data/user/chraab/Output/fast_response/multisample/inherit/')
    
    args = parser.parse_args()

    if '_' in args.name:
        print('Warning: underscores in source name cause LaTeX to fail!')
    tic = datetime.datetime.now()

    args.stop = (Time(args.start, format='iso') + TimeDelta(args.duration, format='jd')).iso
    
    run_trials(args)
    toc = datetime.datetime.now()
    print((toc - tic).total_seconds())

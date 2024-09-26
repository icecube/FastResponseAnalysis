r'''Script to initiate fast reponse
    followups.

    Author: Alex Pizzuto
    Date: 2021
    '''

from fast_response.MultiExternalFollowup import MultiFollowup, GFUFollowup, GrecoFollowup
import argparse
import subprocess
import warnings
import fast_response.web_utils as web_utils
import logging as log
import pyfiglet
import os
import numpy as np

from astropy.time import Time, TimeDelta
import datetime

def calculate_sensitivity(args):

    results = ['sensitivity', 'trials', 'weights']
    outdir = {}
    for _r in results:
        _outdir = os.path.join(args.out, _r)
        if not os.path.exists(_outdir):
            os.makedirs(_outdir)
        outdir[_r] = _outdir
    
    followups = []
    if 'GFU' in args.dataset:
        followups.append(GFUFollowup)
    if 'Greco' in args.dataset:
        followups.append(GrecoFollowup)
    MultiFollowup._followups = followups

    for attr in ['index']:
        setattr(MultiFollowup, f'_{attr}', getattr(args, attr))
    
    
    f = MultiFollowup(
        args.name, args.ra, args.dec, args.start,
        args.stop, extension=args.extension,    
        skipped=args.skip_events, save=False,
    )
    assert f.scramble

    print(str(f))
    print(f.dataset)
    print(args.out)
    print(f.outdir)
    print(f.save_output,
          f.analysisid,
          f.analysispath,
         )
    f.initialize_injector()

    dataset_string = '+'.join(sorted(f.datasets))

    label = f'{f.dec:.03f}_{args.duration}_{f._index}_{dataset_string}'
    outpath = {_r:os.path.join(outdir[_r], f'{label}.npy') for _r in results}
    if os.path.exists(outpath['trials']):
        trials = np.load(outpath['trials'])
        print(f'loading {trials.size} previous trial')
    else:
        trials = None

    sensitivity, trials, weights = f.llh.weighted_sensitivity(0.5, 0.9, f.inj,
                         src_ra = f.ra, src_dec = f.dec,
                         eps=args.eps,
                         n_iter=args.n_iter,
                         n_bckg=args.ntrials,
                         trials = trials,
                         )
    np.save(outpath['trials'], trials)
    np.save(outpath['sensitivity'], sensitivity)
    print(sensitivity)
    

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    log.basicConfig(level=log.ERROR)

    parser = argparse.ArgumentParser(description='Fast Response Analysis')
    parser.add_argument('--name', type=str, default="test inherit sensitivity",
                        help='Name of the source (do not use underscores or LaTeX will be mad!)')
    parser.add_argument('--ra', default=None, type=float,
                        help='Right ascension (in degrees)')
    parser.add_argument('--dec', default=None, type=float,
                        help='Declination (in degrees)')
    parser.add_argument('--duration', default=10, type=float,
                        help='Duration of followup [days]')
    parser.add_argument('--start', type=str, required=False,
                        default=Time(58933.0 - 100, format='mjd').iso,
                        help="Start time of the analysis in ISO format")
    parser.add_argument('--extension', type=float, default=None,
                        help="Source extension in degrees")
    parser.add_argument('--index', type=float, default=None,
                        help="Spectral index to assume in LLH and injected hypothesis")
    parser.add_argument('--skip-events', default=None,
                        type= lambda z:[ tuple(int(y) for y in x.split(':')) for x in z.split(',')],
                        help="Event to exclude from the analyses, eg."
                        "Example --skip-events=127853:67093193")
    parser.add_argument('--ntrials', default=1000, type=int,
                        help="Number of background trials to perform")
    parser.add_argument('--n_iter', default=10, type=int,
                        help="Number of signal trials per per sensitivity iteration")
    parser.add_argument('--eps', default=0.05, type=float,
                        help="Cutoff uncertainty for the sensitivity calculation")
    parser.add_argument('--out', default='/data/user/chraab/Output/fast_response/multisample/inherit/')
    parser.add_argument('--dataset', default=None, type=str,
                        help='Dataset(s) to include, joined by + signs')
    
    args = parser.parse_args()

    if '_' in args.name:
        print('Warning: underscores in source name cause LaTeX to fail!')
    tic = datetime.datetime.now()

    args.stop = (Time(args.start, format='iso') + TimeDelta(args.duration, format='jd')).iso
    
    calculate_sensitivity(args)
    toc = datetime.datetime.now()
    print((toc - tic).total_seconds())

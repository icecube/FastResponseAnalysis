'''Script to calculate sensitivities of analyses defined
in the fast_response framework.

    Author: Christoph Raab
    Date: 2025
    '''

from fast_response.MultiExternalFollowup import MultiFollowup, GFUFollowup, GrecoFollowup, DNNFollowup, DNNOnlineFollowup, DNNIceManFollowup
from fast_response.sensitivity_utils import binomial_error, sensitivity_fit, chi2cdf, find_nearest_idx
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

def grid_scan(f, args, mu_grid, n_iter, TSval):
    sensitivity_trials = {}
    passing_fraction = []
    trials_n = []
    for _mu in mu_grid:
        _trials = f.llh.do_trials(
            n_iter,
            injector=f.inj,
            mean_signal=_mu,
            src_ra = f.ra, src_dec = f.dec,
            )
        sensitivity_trials[_mu] = _trials
        trials_n.append(_trials.size)
        passing_fraction.append(np.count_nonzero(_trials["TS"] > TSval) / _trials.size)
    return sensitivity_trials, trials_n, passing_fraction

def calculate_sensitivity(args):

    results = ['sensitivity_scan', 'trials']
    outdir = {}
    for _r in results:
        _outdir = os.path.join(args.out, _r)
        if 'sensitivity' in _r:
            _outdir = os.path.join(_outdir, f"{args.alpha}_{args.beta}")
        if not os.path.exists(_outdir):
            os.makedirs(_outdir)
        outdir[_r] = _outdir
    
    followups = []
    if 'GFU' in args.dataset:
        followups.append(GFUFollowup)
    if 'Greco' in args.dataset:
        followups.append(GrecoFollowup)
    if 'DNN' in args.dataset:
        followups.append(DNNOnlineFollowup)
    if "IceMan" in args.dataset:
        followups.append(DNNIceManFollowup)
    MultiFollowup._followups = followups

    for attr in ['index', 'fix_index']:
        setattr(MultiFollowup, f'_{attr}', getattr(args, attr))
    MultiFollowup._float_index = not MultiFollowup._fix_index
    
    
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
    if args.eband is None:
        print("NOT USING ENERGY BAND")
        f.initialize_injector()
    else:
        halfwidth = np.sqrt(10**args.ewidth)
        print("USING ENERGIES", (args.eband/halfwidth, args.eband*halfwidth))
        f.initialize_injector(e_range=(args.eband/halfwidth, args.eband*halfwidth))
    

    dataset_string = '+'.join(sorted(f.datasets))
    if args.eband is None:
        label = f'{f.dec:.03f}_{args.duration:.10f}_{f._index:.1f}_{dataset_string}'
    else:
        label = f'{f.dec:.03f}_{args.duration}_{args.eband:.0f}_{dataset_string}'
    outpath = {_r:os.path.join(outdir[_r], f'{label}.npy') for _r in results}
    if os.path.exists(outpath['trials']):
        trials = np.load(outpath['trials'])
        n_bg_loaded = np.count_nonzero(trials["n_inj"] == 0)
        #trials = trials[trials["n_inj"] == 0]
        print(f'loading {n_bg_loaded} previous BG trials')

    else:
        dtype = [("n_inj", int), ("TS", f.llh.dtype_precision)]
        dtype.extend(f.llh.params_dtype()) # no source-specific parameters
        trials = np.empty((0, ), dtype=dtype)
    
    ntrials = args.ntrials - trials[trials["n_inj"] == 0].size
    if ntrials > 0:
        print(f"Drawing an extra {ntrials} BG trials")
        trials = np.append(trials,
                           f.llh.do_trials(n_iter=ntrials,
                                           src_ra = f.ra, src_dec = f.dec,
                                           )[list(trials.dtype.names)])


    

    
    bg_ts = trials["TS"][trials["n_inj"] == 0]
    TSval = np.percentile(bg_ts, 100 * (1 - args.alpha))
    alpha_err = binomial_error(args.alpha, bg_ts.size)
    print(f"TS > {TSval:.2f} in {args.alpha} of BG, err {alpha_err/args.alpha:.2%}")

    
    # initial scan
    print(f"Initial scan {args.mu}")
    mu = np.array(args.mu)
    sensitivity_trials, trials_n, passing_fraction = grid_scan(f, args, mu, args.n_iter//10, TSval)
    # find the cross-over +/- 2
    idx = find_nearest_idx(passing_fraction, args.beta)
    # scan in that region with a finer grid
    _mu = np.linspace(mu[max(0, idx - 2)], mu[min(idx + 3, mu.size-1)], 1 + 4*5)
    print(f"Finer scan over {_mu}")
    _sensitivity_trials, _trials_n, _passing_fraction = grid_scan(f, args, _mu, args.n_iter, TSval)
    # combine the scans
    mu = np.append(mu, _mu)
    sensitivity_trials.update(_sensitivity_trials)
    trials_n += _trials_n
    passing_fraction += _passing_fraction
    
    # sens = sensitivity_fit(
    #     mu,
    #     np.array(passing_fraction),
    #     errs = passing_fraction_unc,
    #     fit_func=chi2cdf,
    #     conf_lev=args.beta,
    # )
    # print(sens)
    flux = f.inj.mu2flux(mu)
    import pandas as pd
    df = pd.DataFrame()
    df['n'] = trials_n
    df["mu"] = mu
    df["fluence"] = 1e9*flux
    df["pf"] = passing_fraction
    df['pf_err'] = binomial_error(df['pf'], df['n'])
    df.sort_values("mu", inplace=True)
    
    print(df)
    simple_sens = np.interp(args.beta, df["pf"], df["mu"])
    print(f"simple interpolation: mu={simple_sens:.3f}, fluence={f.inj.mu2flux(simple_sens)*1e9:.2e} GeV/cm2")
    #df['trials'] = [sensitivity_trials[_mu] for _mu in mu]
    
    
    trials = np.concatenate([trials] + [_v for _v in sensitivity_trials.values()])
    np.save(outpath['trials'], trials)

    result = np.array(
            [(args.alpha, args.beta, TSval, simple_sens, f.inj.mu2flux(simple_sens))],
            dtype=[(k, f.llh.dtype_precision) for k in ["alpha", "TS", "beta", "mu", "flux"]])
    np.save(outpath['sensitivity_scan'], result)
    df.to_parquet(outpath['sensitivity_scan'] + ".parquet")
    

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    log.basicConfig(level=log.ERROR)

    parser = argparse.ArgumentParser(description='Fast Response Analysis sensitivity')
    parser.add_argument('--name', type=str, default="test inherit sensitivity",
                        help='Name of the source (do not use underscores or LaTeX will be mad!)')
    parser.add_argument('--ra', default=None, type=float,
                        help='Right ascension (in degrees)')
    parser.add_argument('--dec', default=None, type=float,
                        help='Declination (in degrees)')
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
    parser.add_argument('--eband', default=None, type=float,
                        help='Determine differential sensitivity within an energy band',
                        )
    parser.add_argument('--ewidth', default=0.5, type=float,
                        help='Energy band width in decades',
                        )
    parser.add_argument('--skip-events', default=None,
                        type= lambda z:[ tuple(int(y) for y in x.split(':')) for x in z.split(',')],
                        help="Event to exclude from the analyses, eg."
                        "Example --skip-events=127853:67093193")
    parser.add_argument('--ntrials', default=100000, type=int,
                        help="Number of background trials to perform")
    parser.add_argument('--n_iter', default=10000, type=int,
                        help="Number of signal trials per per sensitivity iteration")
    parser.add_argument('--out', default='/data/user/chraab/Output/fast_response/multisample/inherit/')
    parser.add_argument('--dataset', default=None, type=str,
                        help='Dataset(s) to include, joined by + signs')
    parser.add_argument("--mu", type=float, nargs="+", default=list(np.arange(0.5, 3.0, 0.2)))
    
    args = parser.parse_args()

    if '_' in args.name:
        print('Warning: underscores in source name cause LaTeX to fail!')
    tic = datetime.datetime.now()

    args.stop = (Time(args.start, format='iso') + TimeDelta(args.duration, format='jd')).iso
    
    calculate_sensitivity(args)
    toc = datetime.datetime.now()
    print((toc - tic).total_seconds())

'''Script to calculate sensitivities of analyses defined
in the fast_response framework.

    Author: Christoph Raab
    Date: 2026
    '''

from fast_response.MultiGWFollowup import DNNOnlineFollowup as MultiFollowup
import argparse
import subprocess
import warnings
import fast_response.web_utils as web_utils
import logging as log
import pyfiglet
import os
import numpy as np
import scipy

from skylab import utils

from astropy.time import Time, TimeDelta
import datetime

# FIXME put this into Skylab instead! Ghastly to have it here
def fit_source(followup, seed, **kwargs):
    llh = followup.llh
    spatial_prior = followup.spatial_prior
    injector = followup.inj
    # fix random seed of injector
    injector.set_rng_seed(seed)

    # this injects _fewer_ events than requested on purpose
    # we want to pretend like we injected the full number
    mean_signal = kwargs.pop("mean_signal", 0)
    if kwargs.pop("poisson", True):
        ni = np.random.poisson(mean_signal)
    else:
        ni = mean_signal
    _, inject = injector.sample(ni, poisson=False)

    # fix random seed of likelihood
    llh.set_rng_seed(seed)

    val = llh.scan(
        0.0,0.0, scramble=True, seed=seed, spatial_prior=spatial_prior,
        time_mask = [followup.duration/2., followup.centertime],
        pixel_scan=[followup.nside, followup._pixel_scan_nsigma],
        inject=inject,
    )
    if val.size == 0:
        ts = 0
    else:
        ts = val['TS_spatial_prior_0'].max()
    return ni, ts


def do_trials(followup, n_iter=10, **kwargs):
    print(".")
    dtype = [("n_inj", int), ("TS", float)]
    
    rng = np.random.RandomState(followup.llh.do_trials_seed)

    results = [fit_source(followup, rng.randint(2**32), **kwargs)
                          for i in range(n_iter)]

    trials = np.empty((n_iter, ), dtype=dtype)

    for i in range(n_iter):
        trials["n_inj"][i] = results[i][0]
        trials["TS"][i] = results[i][1]
    return trials

def weighted_sensitivity(followup, bg_ts, alpha, beta, n_iter=100, trials=[]):
    ts = np.percentile(bg_ts, 100*(1-alpha))
    for _mu in range(1, 5):
        _trials = do_trials(followup, n_iter=n_iter, mean_signal=_mu, poisson=True)
        trials.append(_trials)
    # bg_trials = np.empty((bg_ts.size,), dtype=trials[0].dtype)
    # for i, _ts in enumerate(bg_ts):
    #     bg_trials["n_inj"][i] = 0
    #     bg_trials["TS"][i] = _ts
    # trials.append(bg_trials)
    trials = np.concatenate(trials)
    
    bounds = np.percentile(
    trials["n_inj"][trials["n_inj"] > 0],
    q=[followup.llh._ub_perc, 100. - followup.llh._ub_perc])

    if bounds[0] == 1:
        bounds[0] = np.count_nonzero(trials["n_inj"] == 1) /\
            np.sum(trials["n_inj"] < 2)

    def residual(n):
        return np.log10((utils.poisson_percentile(
            n, trials["n_inj"], trials["TS"], ts)[0] - beta)**2)

    seed = np.argmin([residual(n) for n in np.arange(0., bounds[-1])])

    xmin, fmin, success = scipy.optimize.fmin_l_bfgs_b(
        residual, [seed], bounds=[bounds], approx_grad=True)

    mu = xmin.item()
    b, b_err = utils.poisson_percentile(
                mu, trials["n_inj"], trials["TS"], ts)

    print(
        "Best estimate: mu = {0:.2f} ({1:.2%} +/- {2:.2%})".format(
            mu, b, b_err))
    
    flux = followup.inj.mu2flux(mu)
    print("mu = {0:.2f}, flux = {1:.2e}".format(mu, flux))

    values = [(alpha, ts, beta, mu, flux)]
    result = np.empty(
            (1,),
            dtype=[(k, float) for k in ["alpha", "TS", "beta", "mu", "flux"]])

    result.flat = values
    return result, trials
    

def run_trials(args):
    results = ['prior_sensitivity', 'prior_trials', 'prior_weights']
    outdir = {}
    for _r in results:
        _outdir = os.path.join(args.out, _r)
        if _r == 'prior_sensitivity':
            _outdir = os.path.join(_outdir, f"{args.alpha}_{args.beta}")
        if not os.path.exists(_outdir):
            os.makedirs(_outdir)
        outdir[_r] = _outdir
    
    # followups = []
    # if 'GFU' in args.dataset:
    #     followups.append(GFUFollowup)
    # if 'DNN' in args.dataset:
    #     followups.append(DNNOnlineFollowup)
    # MultiFollowup._followups = followups

    for attr in ['index', 'fix_index']:
        setattr(MultiFollowup, f'_{attr}', getattr(args, attr))
    MultiFollowup._float_index = not MultiFollowup._fix_index
    
    
    f = MultiFollowup(
        args.name, args.skymap, args.start,
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
        f.initialize_injector()
    else:
        halfwidth = np.sqrt(10**args.ewidth)
        f.initialize_injector(e_range=(args.eband/halfwidth, args.eband*halfwidth))

    dataset_string = '+'.join(sorted(f.datasets))
    skymap_string = os.path.split(args.skymap)[1].split('.')[0]
    if args.eband is None:
        label = f'{skymap_string}_{args.duration}_{f._index}_{dataset_string}'
    else:
        label = f'{skymap_string}_{args.duration}_{args.eband:.0f}_{dataset_string}'
    outpath = {_r:os.path.join(outdir[_r], f'{label}.npy') for _r in results}

    if os.path.exists(outpath['prior_trials']):
        trials = [np.load(outpath['prior_trials'])]
        print(f'loading {trials[0].size} previous trial')
    else:
        trials = []
    
    month = Time(args.start).datetime.month
    # FIXME Horrible! do better than hardcoding this
    rate = [
        0.22, 0.22, 0.22, 0.21,
        0.21, 0.21, 0.21, 0.2 ,
        0.22, 0.23, 0.22, 0.23,
        ][month - 1]
    bg_ts = f.load_background_trials(rate=rate, ntrials=args.ntrials)

    sensitivity, trials = weighted_sensitivity(f,
                                               bg_ts, args.alpha, args.beta,
                                               n_iter=args.n_iter,
                                               trials=trials,
                                               )
    

    # FIXME this won't work as it will fit a source at a fixed position
    # sensitivity, trials, weights = f.llh.weighted_sensitivity(args.alpha, args.beta, f.inj,
    #                      eps=args.eps,
    #                      n_iter=args.n_iter,
    #                      n_bckg=args.ntrials,
    #                      trials = trials,
    #                      )
    # print(sensitivity)
    

    np.save(outpath['prior_sensitivity'], sensitivity)
    np.save(outpath['prior_trials'], trials)
    

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
    parser.add_argument('--ntrials', default=10000, type=int,
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
    
    run_trials(args)
    toc = datetime.datetime.now()
    print((toc - tic).total_seconds())

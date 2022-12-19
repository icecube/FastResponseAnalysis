import pycondor
import numpy as np
import argparse
import pwd
import os

parser = argparse.ArgumentParser(
    description='Submit script for precomputing background trials')
parser.add_argument(
    '--tw', type=float, default=1000., #[-0.1, +14]day: 1218240
    help='time window to use (in sec)')
parser.add_argument(
    '--ntrials', type=int, default=30000,
    help='Number of precomputed trials to run (default=30,000)')
args = parser.parse_args()

username = pwd.getpwuid(os.getuid())[0]
if not os.path.exists(f'/scratch/{username}/'):
    os.mkdir(f'/scratch/{username}/')
if not os.path.exists(f'/scratch/{username}/gw/'):
    os.mkdir(f'/scratch/{username}/gw/')
if not os.path.exists(f'/scratch/{username}/gw/condor/'):
    os.mkdir(f'/scratch/{username}/gw/condor')

error = f'/scratch/{username}/gw/condor/error'
output = f'/scratch/{username}/gw/condor/output'
log = f'/scratch/{username}/gw/condor/log'
submit = f'/scratch/{username}/gw/condor/submit'

### Create Dagman to submit jobs to cluster
job = pycondor.Job(
    'gw_precomp_trials',
    '/data/user/jthwaites/FastResponseAnalysis/fast_response/precomputed_background/precompute_ts.py',
    error=error,
    output=output,
    log=log,
    submit=submit,
    getenv=True,
    universe='vanilla',
    verbose=2,
    request_cpus=8,
    request_memory=8000,
    extra_lines=[
        'should_transfer_files = YES',
        'when_to_transfer_output = ON_EXIT']
    )

#scan over bg rates
rate = [6.2, 6.4, 6.6, 6.8, 7.0, 7.2]
for bg in rate:
    for seed in range(int(args.ntrials/100)):
        job.add_arg(f'--deltaT {args.tw} --ntrials 100 --seed {seed} --bkg {bg}')

#dagman = pycondor.Dagman(
#    'gw_precomp_trials',
#    submit=submit, verbose=2)
    
#dagman.add_job(job)
#dagman.build_submit()
job.build_submit()
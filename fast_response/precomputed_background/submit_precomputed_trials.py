import pycondor
import numpy as np
import argparse
import pwd
import os, sys
from astropy.time import Time
import fast_response

default_outdir = os.environ.get('FAST_RESPONSE_OUTPUT') 
if default_outdir is None: default_outdir='./'

parser = argparse.ArgumentParser(
    description='Submit script')
parser.add_argument(
    '--deltaT', type=float, default=1000., #[-0.1, +14]day: 1218240
    help='time window to use (in sec)')
parser.add_argument(
    '--ntrials', type=int, default=30000,
    help='Number of precomputed trials to run, default 30,000')
parser.add_argument(
    '--seed_start', type=int,default=0,
    help='Seed to start with when running trials')
parser.add_argument(
    '--n_per_batch', default=500, type=int,
    help='Number of trials to run in each set (default: 500)')
parser.add_argument(
    '--type', default='alert', type=str,
    help='run alert or GW trials. options are alert (default) or gw')
parser.add_argument(
    '--outdir', default= default_outdir, type=str, 
    help = 'Output directory. default is FAST_RESPONSE_OUTPUT (if set) or ./ if not')
args = parser.parse_args()

print('Output trials to directory: {}'.format(args.outdir))

fra_dir = os.path.dirname(fast_response.__file__)
if args.type=='alert':
    script = os.path.join(fra_dir, 'precomputed_background/precompute_ts.py')
    name = 'alert_precomp_trials'
    deltaT = args.deltaT
elif args.type=='gw':
    script = os.path.join(fra_dir, 'precomputed_background/precompute_ts_gw.py')
    name = 'gw_precomp_trials'
    deltaT = int(args.deltaT)
else:
    print('--type must be alert or gw. Check options and try again')
    sys.exit()

username = pwd.getpwuid(os.getuid())[0]
if not os.path.exists(f'/scratch/{username}/'):
    os.mkdir(f'/scratch/{username}/')
if not os.path.exists(f'/scratch/{username}/fra/'):
    os.mkdir(f'/scratch/{username}/fra/')
if not os.path.exists(f'/scratch/{username}/fra/condor/'):
    os.mkdir(f'/scratch/{username}/fra/condor')

error = f'/scratch/{username}/fra/condor/error'
output = f'/scratch/{username}/fra/condor/output'
log = f'/scratch/{username}/fra/condor/log'
submit = f'/scratch/{username}/fra/condor/submit'

### Create Dagman to submit jobs to cluster 
job = pycondor.Job(
    name,
    script,
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

for bg_rate in [6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2]:
    for seed in range(args.seed_start, int(args.ntrials/args.n_per_batch)+1):
        job.add_arg('--deltaT {} --ntrials {} --seed {} --bkg {} --outdir {}'.format(
            deltaT, args.n_per_batch, seed, bg_rate, args.outdir))

dagman = pycondor.Dagman(
    'fra_precompbg',
    submit=submit, verbose=2)
    
dagman.add_job(job)
dagman.build_submit()
#job.build_submit()


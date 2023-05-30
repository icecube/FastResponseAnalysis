#!/usr/bin/env python

import pycondor
import numpy as np
import argparse
import pwd
import os
from astropy.time import Time
import glob

parser = argparse.ArgumentParser(
    description='Submit script')
parser.add_argument(
    '--tw', type=float, default=1000., #[-0.1, +14]day: 1218240
    help='time window to use (in sec)')
parser.add_argument(
    '--dir', type=str, default='./',
    help='directory where trials are loaded from, saves output to dir+/glob_trials/')
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
    '/data/user/jthwaites/FastResponseAnalysis/fast_response/precomputed_background/glob_precomputed_trials.py',
    error=error,
    output=output,
    log=log,
    submit=submit,
    getenv=True,
    universe='vanilla',
    verbose=2,
    request_cpus=8,
    request_memory=10000,
    extra_lines=[
        'should_transfer_files = YES',
        'when_to_transfer_output = ON_EXIT']
    )

job.add_arg('--deltaT {} --dir {}'.format(args.tw, args.dir))

job.build_submit()


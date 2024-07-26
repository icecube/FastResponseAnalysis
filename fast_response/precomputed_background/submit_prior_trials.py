'''
Script to submit a bunch of trials with random seeds
with a given spatial prior
As needed for bkg trials in realtime
'''

import pycondor
import argparse
import pwd, os
import fast_response

parser = argparse.ArgumentParser(
    description='Submit script')
parser.add_argument(
    '--outdir', type=str, required=True,
    help = 'Output directory (will make and output to a subdirectory in this directory called trials/)')
parser.add_argument(
    '--skymap', type=str, required=True,
    help='Skymap path or link')
parser.add_argument(
    '--alert_mjd', type=float, required=True,
    help='MJD of alert event')
parser.add_argument(
    '--deltaT', type=float, default=172800., #[-0.1, +14]day: 1218240
    help='time window to use (in sec)')
parser.add_argument(
    '--ntrials', type=int, default=10000,
    help='Number of precomputed trials to run, default 30,000')
parser.add_argument(
    '--seed_start', type=int,default=0,
    help='Seed to start with when running trials')
parser.add_argument(
    '--n_per_batch', default=500, type=int,
    help='Number of trials to run in each set (default: 500)')
# parser.add_argument(
#     '--type', default='alert', type=str,
#     help='run alert or GW trials. options are alert (default) or gw')
args = parser.parse_args()

if not os.path.exists(os.path.join(args.outdir, 'trials')):
    os.mkdir(os.path.join(args.outdir, 'trials'))
print('Output trials to directory: {}/trials'.format(args.outdir))

fra_dir = os.path.dirname(fast_response.__file__)
script = os.path.join(fra_dir, 'precomputed_background/alert_bkg_trials_withprior.py')
deltaT = args.deltaT
skymap = args.skymap
alert_mjd = args.alert_mjd

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
    'alert_precomp_trials',
    script,
    error=error,
    output=output,
    log=log,
    submit=submit,
    getenv=True,
    universe='vanilla',
    verbose=2,
    request_cpus=8,#10,
    request_memory=8000,#10000,
    extra_lines=[
        'should_transfer_files = YES',
        'when_to_transfer_output = ON_EXIT']
    )

for trial in range(int(args.ntrials / args.n_per_batch)):
    job.add_arg(f'--deltaT {deltaT} --ntrials {args.n_per_batch} --outdir {args.outdir} --skymap {skymap} --alert_mjd {alert_mjd}')

dagman = pycondor.Dagman(
    'fra_precompbg',
    submit=submit, verbose=2)
    
dagman.add_job(job)
dagman.build_submit()


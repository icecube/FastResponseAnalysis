#!/usr/bin/env python
import itertools
import subprocess
import processing_tools as tools

config = dict(
    script='/data/user/chraab/metaprojects/realtime/FastResponseAnalysis/fast_response/precomputed_background/precompute_ts_multi.py',
    deltaT=1000,
    dataset="DNN",
    ntrials=100,
    bkg=[0.20, 0.21, 0.22, 0.23, 0.24],
    seed=list(range(20001+10000, 20001+100000, 1000)),
    nside=64,
    fix_index='',
    index=2,
    outdir='/data/user/chraab/fast_response/multisample/variable_jitter_new_environment',
)
if 'fix_index' not in config:
    config['outdir'] += '_floating'


with open('jobs.txt', 'w') as f:
    for _config in tools.generate_arg_combinations(config):
        script = _config.pop('script')
        arg_string = tools.build_arg_string(_config)
        f.write(' '.join([script, arg_string, '\n']))

subprocess.run('condor_submit job_alma9.sub', shell=True, check=True)

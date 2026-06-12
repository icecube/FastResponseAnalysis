#!/usr/bin/env python
import itertools
import subprocess
import numpy as np
import processing_tools as tools
dec_grid = [-71.80512766, -58.21166938, -48.59037789, -40.54160187,
       -33.36701297, -26.74368395, -20.48731511, -14.47751219,
        -8.62692656,  -2.86598398,   2.86598398,   8.62692656,
        14.47751219,  20.48731511,  26.74368395,  33.36701297,
        40.54160187,  48.59037789,  58.21166938,  71.80512766]


config = dict(
    script='/data/user/chraab/metaprojects/realtime/FastResponseAnalysis/fast_response/scripts/multisample_sensitivity.py',
    ra=0.,
    #dec=[19.8, -30.],
    dec = dec_grid, index=[2.0, 2.5],
    duration=[1000./3600/24],
    #duration=[1000./3600/24,2,14],
    #dec = [-30, 0., 30], index=[2.5, 2.0],
    #eband=list(10.**np.arange(1,7+1, 0.5)[1::2]),
    dataset=[
            #'Greco',
             #'GFU+Greco',
            #  'GFU',
             #'GFU+DNN',
    #         'DNN',
             'IceMan',
             #'GFU+Greco+DNN',
             ],
    eps=0.01,
    fix_index='',
    start = '2018-10-01',
    #start = '2025-06-01',
    ntrials=100000,
    alpha=1.35e-3, beta=0.5,
    n_iter=300,
    out='/data/user/chraab/Output/fast_response/multisample/variable_jitter_new_environment',
)
if "fix_index" not in config:
    config["out"] += "_floating"

with open('jobs.txt', 'w') as f:
    for _config in tools.generate_arg_combinations(config):
        script = _config.pop('script')
        arg_string = tools.build_arg_string(_config)
        f.write(' '.join([script, arg_string, '\n']))

subprocess.run('condor_submit job_alma9.sub', shell=True, check=True)

#!/usr/bin/env python
import itertools
import subprocess
import numpy as np
import processing_tools as tools
# dec_grid = [-71.80512766, -58.21166938, -48.59037789, -40.54160187,
#        -33.36701297, -26.74368395, -20.48731511, -14.47751219,
#         -8.62692656,  -2.86598398,   2.86598398,   8.62692656,
#         14.47751219,  20.48731511,  26.74368395,  33.36701297,
#         40.54160187,  48.59037789,  58.21166938,  71.80512766]
# dec_grid = np.linspace(-1, 1, )
dec_grid = list(np.arange(-90, 90 + 1, 5)[1:-1])


config = dict(
    script='/data/user/chraab/metaprojects/realtime/FastResponseAnalysis/fast_response/scripts/multisample_sensitivity_scan.py',
    ra=0.,
    dec = dec_grid, index=[2.0],
    duration=[1000./3600/24],
    #eband=list(10.**np.arange(1,7+1, 0.5)[1::2]),
    dataset=[
            #'Greco',
             #'GFU+Greco',
            #  'GFU',
             #'GFU+DNN',
             'DNN',
    #         'IceMan',
             #'GFU+Greco+DNN',
             ],
    start = '2026-06-30',
    ntrials=100000,
    n_iter=10000,
    out='/data/user/chraab/fast_response/multisample/extended_archival',
)
if "fix_index" not in config:
    config["out"] += "_floating"

quantities = [
    dict(alpha=1.35e-3, beta=0.5),
    dict(alpha=0.5, beta=0.9),
]


with open('jobs.txt', 'w') as f:
    for q in quantities:
        for _config in tools.generate_arg_combinations(config|q):
            script = _config.pop('script')
            arg_string = tools.build_arg_string(_config)
            f.write(' '.join([script, arg_string, '\n']))

subprocess.run('condor_submit job_alma9.sub', shell=True, check=True)

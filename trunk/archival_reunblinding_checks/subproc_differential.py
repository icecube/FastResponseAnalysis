import subprocess
import numpy as np

delta_ts = np.logspace(3., 6., 4)
sinDecs = [-0.5, 0.0, 0.5]

for delta_t in delta_ts:
    for sd in sinDecs:
        subprocess.call(['python', './differential_sensitivity.py', '--deltaT={}'.format(delta_t), '--sinDec={}'.format(sd)])

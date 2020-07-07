import subprocess
import numpy as np
import pandas as pd
import argparse

#parser = argparse.ArgumentParser(description='Fast Response Analysis')
#parser.add_argument('--floor', type=float, default=0.2,
#                    help='Angular Uncertainty floor in deg.: 0.5 in grbllh, nu-sources 0.2')
#args = parser.parse_args()

results_df = pd.read_pickle('/data/user/apizzuto/fast_response_skylab/results_dataframe_grbllh.pkl')

for ind in range(len(results_df.index.values)):
    for floor in [0.2]:
        subprocess.call(['python', './reunblind_archival_analyses.py', '--index={}'.format(ind), '--floor={}'.format(floor)])    

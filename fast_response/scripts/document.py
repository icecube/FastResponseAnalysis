r'''Take already run analysis and document it

    Author: Alex Pizzuto
    Date: May, 2020
    '''

import argparse
import subprocess
import pickle
from glob import glob
import fast_response.web_utils as web_utils

parser = argparse.ArgumentParser(description='Document for FRA')
parser.add_argument('--path', type=str,default=None,
                    help='Path to analysis')
args = parser.parse_args()

with open(glob(args.path + '*_results.pickle')[0], 'rb') as f:
    results = pickle.load(f)

subprocess.call(['cp','-r', results['analysispath'],
        '/home/apizzuto/public_html/FastResponse/webpage/output/{}'.format(results['analysisid'])])
web_utils.updateFastResponseWeb(results)

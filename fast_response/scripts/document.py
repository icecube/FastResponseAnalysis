#!/usr/bin/env python

r'''Take already run analysis and document it

    Author: Alex Pizzuto
    Date: May, 2020
    '''

import argparse
import subprocess
import pickle
from glob import glob
import os, pwd
import fast_response.web_utils as web_utils

parser = argparse.ArgumentParser(description='Document for FRA')
parser.add_argument('--path', type=str,default=None,
                    help='Path to analysis')
parser.add_argument('--gw', action = 'store_true', default=False,
                    help='save GW ')
parser.add_argument('--update', action='store_true', default=False,
                    help='Raise to update values (e.g. if more bkg trials are run)')
args = parser.parse_args()

analysis_path = os.path.join(args.path, '')
with open(glob(analysis_path + '*_results.pickle')[0], 'rb') as f:
    results = pickle.load(f)

username = pwd.getpwuid(os.getuid())[0]
if username == 'realtime': username='jthwaites'

if args.gw:
    subprocess.call(['cp','-r', results['analysispath'],
        '/home/{}/public_html/FastResponse/gw-webpage/output'.format(username)])
else: 
    subprocess.call(['cp','-r', results['analysispath'],
        '/home/{}/public_html/FastResponse/webpage/output'.format(username)])

web_utils.updateFastResponseWeb(results, gw=args.gw, update=args.update)

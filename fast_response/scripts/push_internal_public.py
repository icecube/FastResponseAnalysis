#!/usr/bin/env python

import pandas as pd
from fast_response.web_utils import update_internal_public
import argparse, glob, os, sys

parser = argparse.ArgumentParser(description='Fast Response Followup')
parser.add_argument('--alert_circ', type=int, default=None,
                    help='GCN circular number for the initial alert')
parser.add_argument('--alert_notice', default=None,
                    type=lambda z:tuple(int(y) for y in z.split(':')),
                    help ='run_id:event_id of the alert event, to link notice' + \
                        'if no circular was sent. Circular preferred if possible')
parser.add_argument('--fra_circ', type=int, default=None,
                    help='GCN circular number for FRA circular')
parser.add_argument('--alert_name', type=str, default=None,
                    help='Alert identifier. Format - track: IceCube-YYMMDDA, '+ \
                        'cascade: IceCube-Cascade YYMMDDA')
parser.add_argument('--save_location', type=str, default=None,
                    help='Path to saved FRA results (default uses FAST_RESPONSE_OUTPUT env variable)')
args = parser.parse_args()

output = os.environ.get('FAST_RESPONSE_OUTPUT') if args.save_location is None else args.save_location
name = args.alert_name.replace(' ','_')

if args.alert_circ is not None:
    alert_gcn = args.alert_circ
elif args.alert_notice is not None:
    alert_gcn = args.alert_notice
else:
    print('Must either give the Circular number or the Run_id/event_id for the alert event.')
    sys.exit()

results_1000s = glob.glob(os.path.join(output, '*{}_1.0e+03_s/*.pickle'.format(name)))
results_2d = glob.glob(os.path.join(output, '*{}_1.7e+05_s/*.pickle'.format(name)))

if len(results_1000s) != 1 or len(results_2d) != 1:
    print('Results pickle files not found. Check output location and re-run')
    sys.exit()

with open(results_1000s[0], 'rb') as f:
    analysis_1000 = pd.read_pickle(f)
with open(results_2d[0], 'rb') as fi:
    analysis_2d = pd.read_pickle(fi)

update_internal_public(analysis_1000, analysis_2d, 
                       alert_gcn=alert_gcn, fra_circular=args.fra_circ)
print('done.')

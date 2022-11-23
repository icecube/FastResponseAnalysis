#!/usr/bin/env python

r'''Script to run followup in response to 
gravitational wave events from LVC

Author: Raamis Hussain, updated by Jessie Thwaites
April 2022'''

import argparse, subprocess
from astropy.time import Time
import pyfiglet

from fast_response.GWFollowup import GWFollowup
#import fast_response.web_utils as web_utils

parser = argparse.ArgumentParser(description='GW Followup')
parser.add_argument('--skymap', type=str, default=None,
                    help='path to skymap (can be a web link for GraceDB (LVK) or a path')
parser.add_argument('--time', type=float, default=None,
                    help='Time of the GW (mjd)')
parser.add_argument('--name', type=str,
                    default="name of GW event being followed up")
parser.add_argument('--allow_neg_ts', type=bool, default=False,
                    help='bool to allow negative TS values in gw analysis.')
args = parser.parse_args()

#GW message header
message = '*'*80
message += '\n' + str(pyfiglet.figlet_format("GW Followup")) + '\n'
message += '*'*80
print(message)

delta_t = 1000.
gw_time = Time(args.time, format='mjd')
start_time = gw_time - (delta_t / 86400. / 2.)
stop_time = gw_time + (delta_t / 86400. / 2.)
start = start_time.iso
stop = stop_time.iso

name = args.name
name = name.replace('_', ' ')

#import time
#start_counter = time.perf_counter()
f = GWFollowup(name, args.skymap, start, stop)
f._allow_neg = args.allow_neg_ts

#print('\n Time to initialize: %.2f'%(time.perf_counter()-start_counter))
f.unblind_TS()
f.plot_ontime()
#start_bg_dist = time.perf_counter()-start_counter
f.calc_pvalue()
#print('\n time for bg dist and p-val calc: %.2f'%(time.perf_counter()-start_bg_dist))

f.make_dNdE()
f.plot_tsd(allow_neg=f._allow_neg)
f.upper_limit()
f.find_coincident_events()
f.per_event_pvalue()
results = f.save_results()
f.generate_report()

f.write_circular() 

#print('\n Total time: %.2f'%(time.perf_counter()-start_counter))
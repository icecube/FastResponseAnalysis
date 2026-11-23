#!/usr/bin/env python

r'''Script to run followup in response to 
gravitational wave events from LVC

Author: Sam Hori
November 2026'''

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
parser.add_argument('--tw', default = 1000, type=int, #2 week: 1218240
                    help = 'Time window for the analysis, in sec (default = 1000)')
parser.add_argument('--allow_neg_ts', type=bool, default=False,
                    help='bool to allow negative TS values in gw analysis.')
args = parser.parse_args()

#GW message header
message = '*'*80
message += '\n' + str(pyfiglet.figlet_format("GW Followup")) + '\n'
message += '*'*80
print(message)

gw_time = Time(args.time, format='mjd')
delta_t = float(args.tw)

if args.tw==1000:
    start_time = gw_time - (delta_t / 86400. / 2.)
    stop_time = gw_time + (delta_t / 86400. / 2.)
else: #2 week followup
    print('Beginning 2 week NS followup')
    start_time = gw_time - 0.1
    stop_time = gw_time + 14.

start = start_time.iso
stop = stop_time.iso

name = args.name
name = name.replace('_', ' ')

f = GWFollowup(name+" posterior", args.skymap, start, stop)
f.make_posterior_map(fluxToTest=[.25,.5,.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4], 
                    plotting_location=f.analysispath + '/',
                    )
f.generate_posterior_report()

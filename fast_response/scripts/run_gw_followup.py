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
                    help='path to skymap')
parser.add_argument('--time', type=float, default=None,
                    help='Time of the GW (mjd)')
parser.add_argument('--name', type=str,
                    default="GW Followup")
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

f = GWFollowup(name, args.skymap, start, stop)
f.unblind_TS()
f.plot_ontime()
f.calc_pvalue()
f.make_dNdE()
f.plot_tsd()
f.upper_limit()
f.find_coincident_events()
f.per_event_pvalue()
results = f.save_results()
f.generate_report()

f.write_circular() #need to update! to write proper gw options


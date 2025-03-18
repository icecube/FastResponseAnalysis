#!/usr/bin/env python

r'''
Author: Sam Hori
July 2024'''

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
parser.add_argument('--n_inj', default = 0, type=int,
                    help = 'Number to inject')
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

f = GWFollowup(name, args.skymap, start, stop)
t_mask=(f.llh.exp['time']<=f.stop)&(f.llh.exp['time']>=f.start)

print("__________________________________",f.start,f.stop)
print(len(f.llh.exp),f.llh.exp.dtype.names,
        len(f.llh.exp[t_mask]))
print(f.llh.exp['time'])

f._allow_neg = args.allow_neg_ts
f.save_items['merger_time'] = Time(gw_time, format='mjd').iso
f.save_items['skymap_link'] = args.skymap
f.initialize_injector()
injected_events=f.inject_events(args.n_inj)
f.unblind_TS(custom_events=injected_events)
f.make_prior_map(plotting_location=f.analysispath + '/',format=True)

if args.tw>1000:
    f.plot_ontime(label_events = False,custom_events=injected_events)
else:
    f.plot_ontime(custom_events=injected_events)

f.calc_pvalue()
f.make_dNdE()
f.plot_tsd(allow_neg=f._allow_neg)
#if delta_t/86400. > 1.:
#    f.get_best_fit_contour()

f.upper_limit()
f.find_coincident_events(custom_events=injected_events)
f.per_event_pvalue()

results = f.save_results()
nsToTest=[.25,.5,.75,1,1.25,1.5,1.75,2,2.25,2.50,2.75,3,3.25,3.5,3.75,4]
f.make_posterior_map(fluxToTest=nsToTest, 
                     prior_func=lambda ns,gamma,ra,dec: 1, 
                     plotting_location=f.analysispath + '/',
                     custom_events=injected_events)
f.generate_report()
f.generate_posterior_report()

f.write_circular() 


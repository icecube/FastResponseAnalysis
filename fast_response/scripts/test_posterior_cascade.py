#!/usr/bin/env python

r'''Script to run followup to 
a realtime cascade alert

Author: Alex Pizzuto
May 2020'''

import numpy as np
import os, sys, argparse, subprocess
from astropy.time import Time

from fast_response.AlertFollowup import CascadeFollowup
import fast_response.web_utils as web_utils

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--skymap', type=str, default=None,
                    help='path to skymap')
parser.add_argument('--time', type=float, default=None,
                    help='Time of the alert event (mjd)')
parser.add_argument('--gcn_notice_num', default=0, type=int,
                    help="Number of GCN circular. If not set, links to automated notice")
parser.add_argument('--alert_id', default=None,
                    type= lambda z:[ tuple(int(y) for y in x.split(':')) for x in z.split(',')],
                    help="list of events to exclude from this analysis. "
                    "such as HESE events that contributed to the trigger."
                    "Example --alert_id  127853:67093193,128290:6888376")
parser.add_argument('--suffix', type=str, default='A',
                    help="letter to differentiate multiple alerts on the same day (default = A)."
                    "Event name given by IceCube-yymmdd + suffix.")
args = parser.parse_args()

cascade_time = Time(args.time, format='mjd')
year, month, day = cascade_time.iso.split('-')
day = day[:2]
casc_name = 'IceCube-Cascade_{}{}{}{}'.format(year[-2:], month, day, args.suffix)

if 'https://roc.icecube.wisc.edu' in args.skymap:
    print('Downloading skymap from https://roc.icecube.wisc.edu')
    saved_skymaps_path = os.environ.get('FAST_RESPONSE_OUTPUT') + '/../cascade_skymaps/'
    skymap_filename=args.skymap.split('/')[-1]
    if not os.path.isdir(saved_skymaps_path):
        subprocess.call(['mkdir', saved_skymaps_path])
    subprocess.call(['wget', args.skymap,'--no-check-certificate'])
    subprocess.call(['mv',skymap_filename,saved_skymaps_path])
    print('Done.')
    args.skymap=saved_skymaps_path+skymap_filename

all_results = {}
for delta_t in [2.*86400.]: #1000
    start_time = cascade_time - (delta_t / 86400. / 2.)
    stop_time = cascade_time + (delta_t / 86400. / 2.)
    start = start_time.iso
    stop = stop_time.iso

    name = casc_name + ' {:.1e}_s'.format(delta_t)
    name = name.replace('_', ' ')

    run_id = args.alert_id[0][0]
    ev_id = args.alert_id[0][1]

    f = CascadeFollowup(name, args.skymap, start, stop, skipped=args.alert_id)
    t_mask=(f.llh.exp['time']<=f.stop)&(f.llh.exp['time']>=f.start)
    print("__________________________________",f.start,f.stop)
    print(len(f.llh.exp),f.llh.exp.dtype.names,
          len(f.llh.exp[t_mask]))
    print(f.llh.exp['time'])
    print(max(f.llh.exp['time']))
    print(max(f.exp['time']))

    f.unblind_TS()
    f.plot_ontime()
    f.calc_pvalue()
    f.make_dNdE()
    f.plot_tsd()
    f.upper_limit()
    f.find_coincident_events()
    results = f.save_results()
    f.make_prior_map(plotting_location=f.analysispath + '/',format=True)

    f.generate_report()

    nsToTest=[.25,.5,.75,1,1.25,1.5,1.75,2,2.25,2.50,2.75,3,3.25,3.5,3.75,4,4.25,4.50,4.75,5]
    f.make_posterior_map(nsToTest=nsToTest, 
                     prior_func=lambda ns,gamma,ra,dec: 1, 
                     plotting_location=f.analysispath + '/',
                     )


    all_results[delta_t] = results

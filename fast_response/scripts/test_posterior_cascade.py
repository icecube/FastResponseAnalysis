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
parser.add_argument('--prior', type=str, default="constant",
                    help='instructions on the prior to submit')
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

    f = CascadeFollowup(name+"_"+args.prior, args.skymap, start, stop, skipped=args.alert_id)
    if(args.suffix=="high_ene"):
        print("changing energies")
        f.llh.exp["logE"][f.llh.exp["event"]==50721978]=6.5
        f.llh.exp["logE"][f.llh.exp["event"]==43629278]=6.5
    if(args.suffix=="high_bkg"):
        print("adding events")
        print(len(f.exp))
       
        #just add several multiples  to the dataset
        for i in range(6):
            import copy
            temp1=copy.deepcopy(f.exp)
            
            temp1["ra"]=np.random.uniform(low=0.0, high=2*np.pi, size=len(f.exp["ra"]))
            temp1["azimuth"]=np.random.uniform(low=0.0, high=2*np.pi, size=len(f.exp["ra"]))
            
            temp_mask=(temp1['time']<=f.stop)&(temp1['time']>=f.start)
            temp1["time"][temp_mask]=f.stop+10
            f.exp=np.concatenate((f.exp,temp1),axis=0)
            #f.llh.exp["logE"][f.llh.exp["event"]==50721978]=3
            #f.llh.exp["logE"][f.llh.exp["event"]==43629278]=3
        print(len(f.exp))
        f.llh=f.initialize_llh() 
        f.initialize_injector()
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
    #nsToTest=[1e-7,2e-7,3e-7,4e-7,5e-7, 6e-7,7e-7,8e-7]
    
    if (args.prior=="constant"):
        print("Using Constant Prior")
        prior_func=lambda ns,gamma,ra,dec: 1
    elif(args.prior=="jeffries"):
        print("Using Jeffries Prior")
        prior_func=lambda ns,gamma,ra,dec: 1/np.power(ns,1/2)
    else:
        print("Bad Prior")
        exit()

    f.make_posterior_map(fluxToTest=nsToTest, 
                     prior_func=prior_func, 
                     plotting_location=f.analysispath + '/',
                     useFlux=False
                     )

    f.generate_posterior_report()

    all_results[delta_t] = results

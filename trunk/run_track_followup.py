r'''Script to rerun all archival analyses
with the new Fast Response Analysis Framework.
Relies on having a dataframe with the analysis 
information

Author: Alex Pizzuto
May 2020'''

import numpy as np
import os, sys, argparse
from astropy.time import Time
from astropy.time import TimeDelta
import pandas as pd
import subprocess
import pickle

base_path = os.path.join('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/','')
sys.path.append(base_path)

from FastResponseAnalysis import FastResponseAnalysis
import utils

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--skymap', type=str, default=None,
                    help='path to skymap')
parser.add_argument('--time', type=float, default=None,
                    help='Time of the alert event (mjd)')
parser.add_argument('--document', default=False, action='store_true')
parser.add_argument('--gcn_notice_num', default=0, type=int)
args = parser.parse_args()

track_time = Time(args.time, format='mjd')
year, month, day = track_time.iso.split('-')
day = day[:2]
track_name = 'IceCube-{}{}{}A'.format(year[-2:], month, day)

all_results = {}
for delta_t in [1000., 2.*86400.]:
    start_time = track_time - (delta_t / 86400. / 2.)
    stop_time = track_time + (delta_t / 86400. / 2.)
    start = start_time.iso
    stop = stop_time.iso

    name = track_name + ' {:.1e}_s'.format(delta_t)

    source = {}
    source['Name'] = name.replace('_', ' ')
    source['alert_event'] = True
    source['smear'] = True
    f = FastResponseAnalysis(args.skymap, start, stop, **source)
    f.unblind_TS()
    f.plot_ontime()
    f.calc_pvalue()
    f.make_dNdE()
    f.plot_tsd()
    f.upper_limit()
    results = f.save_results()
    f.generate_report()
    if args.document:
        subprocess.call(['cp','-r',results['analysispath'],
        '/home/apizzuto/public_html/FastResponse/webpage/output/{}'.format(results['analysisid'])])
        utils.updateFastResponseWeb(results)
    all_results[delta_t] = results

all_results[1000.]['gcn_num'] = args.gcn_notice_num
utils.write_alert_gcn(all_results)


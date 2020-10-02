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

from fast_response.FastResponseAnalysis import FastResponseAnalysis
import fast_response.utils

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--index', type=int,default=None,
                    help='Analysis number in dataframe')
parser.add_argument('--deltat', type=float, default=None,
                    help='Time window for analysis (1000 s or 172800 s)')
args = parser.parse_args()

output_dir = "/data/user/apizzuto/fast_response_skylab/fast-response/fast_response/cascades_results/"
with open('/home/tgregoire/hese_cascades/archival_data.pkl', 'r') as f:
    cascades_archival = pickle.load(f)

src_index = args.index
cascade_time = Time(cascades_archival['mjd'][src_index], format='mjd')
source_name = 'IceCube-Cascade_{}_{}'.format(int(cascades_archival['run'][src_index]), int(cascades_archival['event'][src_index]))
start_time = cascade_time - (args.deltat / 86400. / 2.)
stop_time = cascade_time + (args.deltat / 86400. / 2.)
start = start_time.iso
stop = stop_time.iso

loc = output_dir + 'skymaps/' + source_name + '.hp'
name = source_name + ' {:.1e}_s'.format(args.deltat)

source = {}
source['Name'] = name.replace('/', ' ')
source['Name'] = source['Name'].replace('_', ' ')
source['Name'] = source['Name'] + ' test'
source['floor'] = np.radians(0.2)
source['output_dir'] = output_dir
source['alert_event'] = True
f = FastResponseAnalysis(loc, start, stop, **source)
f.unblind_TS()
f.plot_ontime()
f.calc_pvalue()
f.make_dNdE()
f.plot_tsd()
f.upper_limit()
results = f.save_results()
f.generate_report()
#subprocess.call(['cp','-r',results['analysispath'],
#        '/home/apizzuto/public_html/FastResponse/webpage/output/{}'.format(results['analysisid'])])
#if args.floor == 0.2:
#    utils.updateFastResponseWeb(results)

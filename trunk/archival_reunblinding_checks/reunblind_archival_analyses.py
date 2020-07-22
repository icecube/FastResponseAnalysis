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

base_path = os.path.join('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/','')
sys.path.append(base_path)

from FastResponseAnalysis import FastResponseAnalysis
import utils

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--index', type=int,default=None,
                    help='Analysis number in dataframe')
parser.add_argument('--floor', type=float, default=0.2,
                    help='Angular Uncertainty floor in deg.: 0.5 in grbllh, nu-sources 0.2')
parser.add_argument('--ntrials', type=int, default=1000,
                    help="Number of background trials to run, default 1000")
args = parser.parse_args()

output_dir = "/data/user/apizzuto/fast_response_skylab/fast-response/trunk/archival_unblinding_results/"
results_df = pd.read_pickle('/data/user/apizzuto/fast_response_skylab/results_dataframe_grbllh.pkl')

src_index = args.index
source_name = results_df.index.values[src_index]

iso_start = results_df['Beginning Observation Time (UTC)'][src_index]
iso_stop = iso_start + results_df['Duration'][src_index]

iso_start = iso_start.strftime('%Y-%m-%d %H:%M:%S.%f')
iso_stop = iso_stop.strftime('%Y-%m-%d %H:%M:%S.%f')
#start = dataset.iso2mjd(iso_start)
#stop = dataset.iso2mjd(iso_stop)

dec = results_df['Source Dec (degrees)'][src_index]
ra = results_df['Source RA (degrees)'][src_index]
loc = "{}, {}".format(ra, dec)

source = {}
source['Name'] = source_name.replace('/', ' ')
source['Name'] = source['Name'].replace('_', ' ')
source['floor'] = np.radians(args.floor)
if args.floor != 0.2:
    source['Name'] = source['Name'] + ' floor {:.1f}'.format(args.floor) 

if args.floor == 0.5:
    output_dir = output_dir + 'floor_0.5/'
source['output_dir'] = output_dir

if results_df['Skipped Events'][src_index] is not None:
    source['Skipped Events'] = [tuple(int(x) for x in results_df['Skipped Events'][src_index].split(':'))]

if results_df['Source extension (degrees)'][src_index] != 0.0:
    source['Extension'] = results_df['Source extension (degrees)'][src_index]

#print(source)
#sys.exit()
f = FastResponseAnalysis(loc, iso_start, iso_stop, **source)
f.unblind_TS()
f.plot_ontime()
f.ns_scan()
f.calc_pvalue(ntrials = int(args.ntrials), run_anyway=True)
f.make_dNdE()
f.plot_tsd()
f.upper_limit()
results = f.save_results()
f.generate_report()
subprocess.call(['cp','-r',results['analysispath'],
        '/home/apizzuto/public_html/FastResponse/webpage/output/{}'.format(results['analysisid'])])
if args.floor == 0.2:
    utils.updateFastResponseWeb(results)

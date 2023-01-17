#!/usr/bin/env python

r'''Script to run track followup in response to 
high-energy neutrino alert events

Author: Alex Pizzuto
May 2020'''

import os, argparse
from astropy.time import Time

from fast_response.AlertFollowup import TrackFollowup
import fast_response.web_utils as web_utils

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--skymap', type=str, default=None,
                    help='path to skymap')
parser.add_argument('--time', type=float, default=None,
                    help='Time of the alert event (mjd)')
parser.add_argument('--gcn_notice_num', default=0, type=int,
                    help="GCN Circular NUMBER for updated circular (not Bacodine/Notice)")
parser.add_argument('--alert_id', default=None,
                    type= lambda z:[ tuple(int(y) for y in x.split(':')) for x in z.split(',')],
                    help="list of events to exclude from this analysis. "
                    "such as HESE events that contributed to the trigger."
                    "Example --alert_id  127853:67093193,128290:6888376")
parser.add_argument('--suffix', type=str, default='A',
                    help="letter to differentiate multiple alerts on the same day (default = A)."
                    "Event name given by IceCube-yymmdd + suffix.")
args = parser.parse_args()

track_time = Time(args.time, format='mjd')
year, month, day = track_time.iso.split('-')
day = day[:2]
track_name = 'IceCube-{}{}{}{}'.format(year[-2:], month, day, args.suffix)

all_results = {}
for delta_t in [1000., 2.*86400.]:
    start_time = track_time - (delta_t / 86400. / 2.)
    stop_time = track_time + (delta_t / 86400. / 2.)
    start = start_time.iso
    stop = stop_time.iso

    name = track_name + ' {:.1e}_s'.format(delta_t)
    name = name.replace('_', ' ')

    run_id = args.alert_id[0][0]
    ev_id = args.alert_id[0][1]

    # look for contour files
    base_skymap_path = '/home/followup/output_plots/'
    contour_fs = [base_skymap_path \
            + f'run{run_id:08d}.evt{ev_id:012d}.HESE.contour_{containment}.txt' \
            for containment in ['50', '90']]
    contour_fs = [f for f in contour_fs if os.path.exists(f)]
    if len(contour_fs) == 0:
        contour_fs = None

    f = TrackFollowup(name, args.skymap, start, stop, skipped=args.alert_id)

    f.unblind_TS()
    f.plot_ontime(contour_files=contour_fs)
    f.calc_pvalue()
    f.make_dNdE()
    f.plot_tsd()
    f.upper_limit()
    f.find_coincident_events()
    results = f.save_results()
    f.generate_report()
    all_results[delta_t] = results

all_results[1000.]['gcn_num'] = args.gcn_notice_num

# Write circular to the output directory of the 2 day analysis
f.write_circular(all_results)


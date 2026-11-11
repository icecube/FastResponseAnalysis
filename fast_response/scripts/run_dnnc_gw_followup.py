#!/usr/bin/env python

"""Script to run followup on gravitational wave events
using DNN Cascades.

Author: Christoph Raab, adapted from run_gw_followup.py.
Date: 2026
"""

import argparse
from astropy.time import Time, TimeDelta
import pyfiglet
import datetime

from fast_response.MultiGWFollowup import DNNOnlineFollowup
import logging
import os
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description='GW Followup')
parser.add_argument('--skymap', type=str, default=None,
                    help='path to skymap (can be a web link for GraceDB (LVK) or a path')
parser.add_argument('--time', type=float, default=None,
                    help='Time of the GW (mjd)')
parser.add_argument('--name', type=str,
                    default="name of GW event being followed up")
parser.add_argument('--tw', default = 1000, type=int, choices={1000},
                    help = 'Time window for the analysis, in sec (default = 1000)')
parser.add_argument('--allow_neg_ts', type=bool, default=False,
                    help='bool to allow negative TS values in gw analysis.')
parser.add_argument("--seed", type=int, default=1,
                    help="RNG seed in case of a scrambled test unblinding")
args = parser.parse_args()

def run_followup(args):
    #GW message header
    message = '*'*80
    message += '\n' + str(pyfiglet.figlet_format("DNNC GW Followup")) + '\n'
    message += '*'*80
    print(message)

    # determine symmetric time window
    gw_time = Time(args.time, format='mjd')
    delta_t = TimeDelta(args.tw, format="sec")
    start_time = gw_time - delta_t / 2
    stop_time = gw_time + delta_t / 2
    start = start_time.iso
    stop = stop_time.iso

    # replace underscores for plots
    name = args.name
    name = name.replace('_', ' ')

    f = DNNOnlineFollowup(name, args.skymap, start, stop, seed=args.seed)
    f._allow_neg = args.allow_neg_ts
    f.save_items['merger_time'] = Time(gw_time, format='mjd').iso
    f.save_items['skymap_link'] = args.skymap

    f.unblind_TS()
    f.plot_ontime(plot_zoom=False)
    # load precomputed BG
    llh = f.llh._samples[0] # have 1 sample
    # use the estimate provided by the temporal model
    # from the background window(s): (temporal model has estimated it from +/- 60 days)
    rate = round(1000 * llh.nbackground / (llh.on_livetime * 86400), 2)
    f.load_background_trials(rate=rate)

    # use the above to calculate p-value
    f.calc_pvalue()
    
    # make plots for internal report
    f.make_dNdE()
    f.plot_tsd(allow_neg=f._allow_neg)

    # calculations for GCN notice
    f.upper_limit()
    f.find_coincident_events()
    f.per_event_pvalue()
    
    # save output
    f.save_results()
    #f.generate_report()
    #f.write_circular()

if __name__ == "__main__":
    args = parser.parse_args()
    run_followup(args)
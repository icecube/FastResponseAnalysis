#!/usr/bin/env python

import subprocess
import glob
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--number_of_mocks', default=5, type=int,
                        help='Limit number of mocks being run.')

args = parser.parse_args()

MockGallery = sorted(glob.glob('/data/user/jthwaites/o4-mocks/*/*.xml'))[::-1]

for file in MockGallery[0:args.number_of_mocks]:
    subprocess.call(['/data/user/mromfoe/software/fastresponse/fast_response/listeners/gw_gcn_listener.py', 
                     '--test_path={}'.format(file),
                     '--log_path={}'.format('/home/mromfoe/public_html/logfile.log')])
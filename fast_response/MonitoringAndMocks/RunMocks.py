#!/usr/bin/env python

import subprocess
import glob
import argparse
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('--number_of_mocks', default=5, type=int,
                        help='Limit number of mocks being run.')

args = parser.parse_args()

MockGallery = sorted(glob.glob('/data/user/jthwaites/o4-mocks/*/*.xml'))[::-1]

pieces = [string.split('-') for string in MockGallery]

Puzzle_Box = {}

for prefix in np.unique([pieces[i][1] for i in range(len(pieces))]):
        Puzzle_Box[prefix]=[[],[],[],[]]
for bit in pieces:
        Puzzle_Box[bit[1]][0].append(int(bit[2]))
        Puzzle_Box[bit[1]][1].append(bit[3])
        Puzzle_Box[bit[1]][2].append(bit[0])
        Puzzle_Box[bit[1]][3].append(bit[4]+'-'+bit[5])

answer = []
for key in Puzzle_Box.keys():
        high_file=max(Puzzle_Box[key][0])
        for i in range(len(Puzzle_Box[key][0])):
                if Puzzle_Box[key][0][i]==high_file:
                        answer.append(Puzzle_Box[key][2][i]+'-'+key+'-'+str(high_file)+
                        '-'+Puzzle_Box[key][1][i]+'-'+Puzzle_Box[key][3][i])

for file in answer[0:args.number_of_mocks]:
    subprocess.call(['/data/user/mromfoe/software/fastresponse/fast_response/listeners/gw_gcn_listener.py', 
                     '--test_path={}'.format(file),
                     '--log_path={}'.format('/data/user/mromfoe/software/fastresponse/fast_response/listeners/mock_logfile.log')])

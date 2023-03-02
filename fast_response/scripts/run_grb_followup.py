#!/usr/bin/env python

r'''Script to run followup in response to 
gravitational wave events from LVC

Author: Raamis Hussain, updated by Jessie Thwaites
April 2022'''

import argparse, subprocess
from astropy.time import Time
import pyfiglet
import os
import sqlite3
import pandas
import sys
import numpy as np
import healpy as hp
from mdutils.mdutils import MdUtils
from mdutils import Html
import datetime

import skylab
from skylab.datasets import Datasets
from icecube import icetray
from os.path import expanduser

from fast_response.GRBFollowup import GRBFollowup
#import fast_response.web_utils as web_utils

parser = argparse.ArgumentParser(description='GW Followup')
#parser.add_argument('--skymap', type=str, default=None,
#                    help='path to skymap (can be a web link for GraceDB (LVK) or a path')
#parser.add_argument('--time', type=float, default=None,
#                    help='Time of the GRB (mjd)')
parser.add_argument('--name', type=str,
                    default="name of GRB event being followed up")
parser.add_argument('--allow_neg_ts', type=bool, default=False,
                    help='bool to allow negative TS values in gw analysis.')
args = parser.parse_args()

# Load the database with the sqlite3 module
#if os.path.exists('GRBweb2.sqlite'):
#	os.remove('GRBweb2.sqlite')
#os.system("wget https://icecube.wisc.edu/~grbweb_public/GRBweb2.sqlite")
db = sqlite3.connect('GRBweb2.sqlite')
Summary_table = pandas.read_sql_query("SELECT * from Summary", db)
print(Summary_table.keys())

test_grb = Summary_table[Summary_table.GRB_name==args.name]
print('GRB name', pandas.Series.tolist(test_grb.GRB_name))
#print('Grb datatable t90 start',test_grb.T90_start)
#print('GRB datatable t0 start', test_grb.T0)
test_grb = test_grb.to_records()

print(test_grb)

GBM_truth_1 = str(test_grb.GBM_located)[1]
try:
	GBM_truth = float(GBM_truth_1)
except:
	GBM_truth = 0
print(GBM_truth)
if GBM_truth == 1:
	print('This is the gbm truth', GBM_truth)
	GBM_name = test_grb.GRB_name_Fermi
	GBM_number = str(GBM_name)[5:-2]
	print(GBM_name,GBM_number)
	skymap = f'/data/user/rpmurphy/realtime/healpix/glg_healpix_all_bn{GBM_number}_v00.fit'
	#print(skymap)

elif os.path.exists(f'/data/user/rpmurphy/healpix/{args.name}.fits'):
	skymap = f'/data/user/rpmurphy/healpix/{args.name}.fits'
	
else:
	dec_1 = str(test_grb.decl)[1:-1]
	print('this is the dec', dec_1)
	dec = np.radians(float(dec_1))
	print('this is the ra', test_grb.ra)
	ra_1 = str(test_grb.ra)[1:-1]
	ra = np.radians(float(ra_1))
	source_1 = str(test_grb.pos_error)[1:-1]
	source_uncertainty = np.radians(float(test_grb.pos_error))
	
	theta = np.pi/2. - dec                                  # Convert to healpy stuff
	phi = ra
	nside = 128                                             # Pass as argument if needed
	probs = np.zeros(hp.pixelfunc.nside2npix(nside))        # Array of zeros with length from nside
	pixel = hp.pixelfunc.ang2pix(nside, theta, phi)         # Find the pixel that matches the location of the GRB
	vec = hp.pixelfunc.pix2vec(nside, pixel)                # Find the vector to that pixel

        # nan may not matter anymore
	if np.isnan(source_uncertainty) == True or source_uncertainty < np.radians(1.0):
		source_uncertainty = np.radians(1.0)            # Scan minimum of 1 degree around source

	pixels = hp.query_disc(nside, vec, source_uncertainty)  # Find all the pixels within appropriate degrees of that pixel
	probs[pixels] = 1.

	skymap = f'/data/user/rpmurphy/healpix/{args.name}.fits'

	#if not os.path.exists:
	hp.fitsfunc.write_map(skymap,probs)	

#else:
#	print('Skymap error')

t90_start = test_grb.mjd + test_grb.T90_start/86400 - test_grb.T0/86400

t0_start = test_grb.mjd + test_grb.T0/86400

t0 = Time(t0_start, format='mjd')

print(t0.mjd, t0.isot, t0.iso)

print('this is the t0 start', test_grb.T0)

print('half t90', test_grb.T90/86400/2)

t90_center = t90_start + test_grb.T90/86400/2

print('this is the t90 center', t90_center)

'''
grb_name = test_grb.GRB_name
test_grb_str = str(test_grb.T90_start)
grb_date = '20{}-{}-{}T{}:{}:{}'.format(args.name[3:5],args.name[5:7], args.name[7:9],
		test_grb_str[1:-9],test_grb_str[-9:-7],test_grb_str[-7:-1])
print('grb dat utc',grb_date)
grb_date_mjd = Time(grb_date, format='isot', scale='utc').mjd  #added 221.23 for the T90 timewindow  #mjd for random time #test_grb.mjd
print(np.dtype(grb_date_mjd))
t90_start = grb_date_mjd + test_grb.T90_start/86400 - test_grb.T0
print('mjd t90 start?',t90_start)
GBM= test_grb.GBM_located
#source_uncertainty = test_grb.pos_error
#grb_name, t100_start, t100_stop, ra, dec, gbm_filename, GBM, source_uncertainty = get_grb_info(grb_index,grb_dict=grb_dict)
#grb_name, t100_start, t100_stop, ra, dec, gbm_filename, GBM, source_uncertainty = get_grb_info(grb_index,grb_dict=grb_dict)
'''
#sys.exit()

#GW message header
message = '*'*80
message += '\n' + str(pyfiglet.figlet_format("GRB Followup")) + '\n'
message += '*'*80
print(message)

tw_list = [10] #, 25, 50, 100]#, 250., 500., 1000., 5000., 172800., 1296000.]
ts_list = []
ns_list = []
p_vals = []
sigmas = []
now = Time(datetime.datetime.now()+datetime.timedelta(hours=6),scale='utc',precision=0)
dataset= Datasets['GFUOnline_v001p02']

mdFile = MdUtils(file_name='{}_report'.format(args.name), title='IceCube Fast-Response Report for {}'.format(args.name))
mdFile.new_header(level=1, title='Overview')
mdFile.new_header(level=2, title = 'On-Time Data')
mdFile.new_table(columns=2, rows=4, text=["","","Access Method", "database", "Stream", "neutrino", "Query Time", str(now.strftime('%Y-%m-%d %H:%M:%S'))])
#Add start and stop time to On-time data?
mdFile.new_header(level=2, title='Skylab Analysis Information')
mdFile.new_table(columns = 2, rows=6, text =["","","Skylab Verson", skylab.__version__, "IceTray Path", str(icetray.__path__).replace('_', '\_'), "Created by", expanduser('~')[6:], "Dataset Used", str(dataset.subdir).replace('_',' '), "Dataset Details", str(dataset.name)])


for tw in tw_list:
	delta_t = tw
	print('Running the',tw,'s time window')
	#gw_time = Time(args.time, format='mjd')
	#start_time = gw_time - (delta_t / 86400. / 2.)
	#stop_time = gw_time + (delta_t / 86400. / 2.)

	if tw == 1296000.:
		t90_center_time = Time(t90_center, format='mjd')
		start_time = t90_center_time - 1
		stop_time = t90_center_time + 14
		print(start_time, stop_time)
	else:
		t90_center_time = Time(t90_center, format='mjd')
		start_time =t90_center_time - (delta_t /86400. / 2.)
		stop_time = t90_center_time + (delta_t / 86400. / 2.)
		print(start_time, stop_time)
		
	start_ndarray = start_time.iso
	#start = start_ndarray.tostring()
	start = str(start_ndarray)[2:-2]
	stop_ndarray = stop_time.iso
	#stop = stop_ndarray.tostring()
	stop = str(stop_ndarray)[2:-2]
	name = args.name
	name = name.replace('_', ' ')
	tw = tw
	print(name, skymap, start, stop, tw)
	
	if tw == 172800:
		tw_1 = '2 day'
	elif tw == 1296000:
		tw_1 = '15 day'
	else:
		tw_1 = str(tw) + 'second'
	mdFile.new_header(level=1, title='{} Time Window'.format(tw_1))
	mdFile.new_header(level=2, title='Source Information')
	mdFile.new_table(columns = 2, rows=5, text= ["Source Name", "{}".format(args.name), "Skymap","{}".format(skymap),"Trigger Time","{} (MJD={})".format(t0.iso, t0.mjd), "Start Time", "{} (Trigger-{}s)".format(Time(t0.mjd+tw/2, format= 'mjd').iso,tw/2), "Stop Time", "{} (Trigger+{}s)".format(Time(t0.mjd-tw/2, format='mjd').iso,tw/2) ] )	

	f = GRBFollowup(name, skymap, start, stop, tw)
	f._allow_neg = args.allow_neg_ts
	f._dataset = 'GFUOnline_v001p02'
	f._fix_index = True
	f._float_index = False
	f._season_names = ['IC86, 2017', 'IC86, 2018', 'IC86, 2019']
	
	print('unblind TS')
	ts, ns = f.unblind_TS()
	ts_list.append(ts)
	ns_list.append(ns)
	print('plot ontime')
	f.plot_ontime()
	print('calc_pvalue')
	#f.calc_pvalue()
	p, sigma = f.calc_pvalue()
	try:
		p_vals.append(p)
		sigmas.append(sigma)
	except:
		p_vals.append(1.0)
		sigmas.append(0)
	print(p_vals)
	print('make dNdE')
	f.make_dNdE()
	print('plot tsd')
	f.plot_tsd(allow_neg=f._allow_neg)
	print('upper limit')
	f.upper_limit()
	print('find coincident events')
	f.find_coincident_events()
	print('per event pvalue')
	f.per_event_pvalue()
	results = f.save_results()
	f.generate_report()

'''
f = GWFollowup(name, args.skymap, start, stop)
f.generate_report()
f.write_circular()''' 

mdFile.create_md_file()


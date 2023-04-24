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
import markdown
from markdown.extensions.tables import TableExtension

import skylab
from skylab.datasets import Datasets
from icecube import icetray
from os.path import expanduser

#from fast_response.make_ontime_plots   import make_rate_plots
import icecube.realtime_tools.live

import json

use_urllib2 = True
#import urllib2, urllib

try:
    import urllib2, urllib

except:
    use_urllib2 = False
    import urllib.request, urllib.parse, urllib.error


from fast_response.GRBFollowup import GRBFollowup
#import fast_response.web_utils as web_utils

import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 
#mpl.rcParams.update(mpl.rcParamsDefault)
#mpl.rcParams['text.usetex'] = True

print(sys.path, os.environ)

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


######################
#                    #
#      Functions     #
#                    #
######################


def query_db_runs( time_window):
	# first query the rundata base with the times of the analysis
	# plus 8 hours on either side
	run_url = 'https://live.icecube.wisc.edu/run_info/'
	run_query = {'user':'icecube', 'pass':'skua',
			'start':(Time(start_trigger, format='mjd').iso),
			'stop': (Time(stop_trigger, format='mjd').iso)}
	now = Time(datetime.datetime.now()+datetime.timedelta(hours=6),scale='utc',precision=0)
	if use_urllib2:
		run_table = json.loads(urllib2.urlopen(
		urllib2.Request(run_url, urllib.urlencode(run_query)),
		timeout=500).read())
	else:
		req_data = urllib.parse.urlencode(run_query).encode("utf-8")
		req = urllib.request.Request(run_url)
		with urllib.request.urlopen(req, data=req_data, timeout=500) as fi:
			run_table = json.load(fi)

	#add duration to run table
	for run in run_table:
		runstart = Time(run['start'],format='iso',scale='utc')
		if run['stop'] is None:
			runstop  = now
			run['stop']=runstop.iso
		else:
			runstop  = Time(run['stop'],scale='utc',precision=0)

		dt=int((runstop-runstart).sec)
		hours = int(dt)//3600
		minutes = int(dt)//60-hours*60
		seconds = dt%60

		livetime = (min(time_window[1],runstop)-
			max(time_window[0],runstart)).sec
		run['livetime'] = livetime if livetime > 0 else 0
		if run['livetime'] > 86400. * 2:	
			run['livetime'] = 0
			run['stop'] = run['start']
			dt=0
			hours = 0
			minutes = 0
			seconds = 0

		run['duration_sec'] = dt
		run['duration'] = "{:01d}:{:02d}:{:02d}".format(hours,minutes,seconds)
		run['OK'] = ("OK" if run['status'] in ['SUCCESS']
			and run['lightmode'] in ['dark']
			and run['filter_mode'] in ['PhysicsFiltering']
			and run['run_mode'] in ['PhysicsTrig']
			else "NotOK"
			)

	return run_table

def ontime_table(query_dict):
	newdict=[]
	for event in query_dict:
		newevent = event['value']['data']
		for key,val in list(newevent['reco']['splinempe'].items()):
			newevent['splinempe_'+key]=val
		if Time(newevent['eventtime'],scale='utc',format='iso'):
			newevent['muex']= newevent['reco']['energy']['mpe_muex']
		del newevent['reco']
		newdict.append(newevent)
	events = pandas.DataFrame(newdict)

	if len(events):
		t=Time(list(events['eventtime']),scale='utc',format='iso')
		events['t']=t.datatime
		events['mjd']=t.mjd
	return events


#report_path ='/home/rpmurphy/public_html/{}'.format(args.name)
report_path = '/data/user/rpmurphy/realtime/fast_response/FastResponseAnalysis/fast_response/scripts/{}'.format(args.name)

os.environ['PATH'] += '/usr/bin/'

if os.path.isdir(report_path):
	print("Directory {} already exists. Deleting it ...".format(report_path))
	subprocess.call(['rm', '-r', report_path])
	subprocess.call(['mkdir', report_path])
	# sys.exit()
else:
	subprocess.call(['mkdir', report_path])

# Load the database with the sqlite3 module
#if os.path.exists('GRBweb2.sqlite'):
#	os.remove('GRBweb2.sqlite')
#os.system("wget https://icecube.wisc.edu/~grbweb_public/GRBweb2.sqlite")
db = sqlite3.connect('GRBweb2.sqlite')
Summary_table = pandas.read_sql_query("SELECT * from Summary", db)
print(Summary_table.keys())

test_grb = Summary_table[Summary_table.GRB_name==args.name]
if test_grb.empty:
	grb_name = str(args.name + '*')
	test_grb = Summary_table[Summary_table.GRB_name==grb_name]	

print('GRB name', pandas.Series.tolist(test_grb.GRB_name))
#print('Grb datatable t90 start',test_grb.T90_start)
#print('GRB datatable t0 start', test_grb.T0)
test_grb = test_grb.to_records()

print(test_grb)

print(mpl.rcParams['text.usetex'], '2')

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
	print('healpix map exists')
	skymap = f'/data/user/rpmurphy/healpix/{args.name}.fits'
	print(skymap)	
else:
	print('generating healpix map')
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

tw_list = [10, 25, 50, 100, 250., 500., 1000., 5000., 172800., 1296000.]
ts_list = []
ns_list = []
p_vals = []
sigmas = []
now = Time(datetime.datetime.now()+datetime.timedelta(hours=6),scale='utc',precision=0)
dataset= Datasets['GFUOnline_v001p02']

mdFile = MdUtils(file_name='{}_report'.format(args.name), title='IceCube Fast-Response Report for {}'.format(args.name))
mdFile.new_line('For Internal Use Only')
mdFile.new_header(level=1, title='Overview')
mdFile.new_header(level=2, title = 'On-Time Data')
mdFile.new_table(columns=2, rows=4, text=["","","Access Method", "database", "Stream", "neutrino", "Query Time", str(now.strftime('%Y-%m-%d %H:%M:%S'))])
#Add start and stop time to On-time data?
mdFile.new_header(level=2, title='Skylab Analysis Information')
mdFile.new_table(columns = 2, rows=6, text =["","","Skylab Verson", skylab.__version__, "IceTray Path", str(icetray.__path__).replace('_', '\_'), "Created by", expanduser('~')[6:], "Dataset Used", str(dataset.subdir).replace('_',' '), "Dataset Details", str(dataset.name)])

print(mpl.rcParams['text.usetex'], '1')

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
	start2stop = (start_time, stop_time)
	print(name, skymap, start, stop, tw)
	
	if tw == 172800:
		tw_1 = '2 day'
	elif tw == 1296000:
		tw_1 = '15 day'
	else:
		tw_1 = str(tw) + 'second'
	
	#print(t0.mjd-tw/2, t0.mjd+tw/2, Time((t0.mjd-tw/2), format='mjd').iso, Time(t0.mjd+tw/2, format='mjd').iso)	
	
	start_trigger = t0.mjd-tw/2/86400 ###CLARIFY TIRGGER TIME T0 or T90 CENTEER???		
	stop_trigger = t0.mjd+tw/2/86400
	print(start_trigger, stop_trigger, Time(start_trigger, format='mjd').iso, Time(stop_trigger, format='mjd').iso)

	mdFile.new_header(level=1, title='{} Time Window'.format(tw_1))
	mdFile.new_header(level=2, title='Source Information')
	mdFile.new_table(columns = 2, rows=5, text= ["Source Name", "{}".format(args.name), "Skymap","{}".format(skymap),"Trigger Time","{} (MJD={})".format(t0.iso, t0.mjd), "Start Time", "{} (Trigger-{}s)".format(Time(start_trigger, format= 'mjd').iso, tw/2), "Stop Time", "{} (Trigger+{}s)".format(Time(stop_trigger, format='mjd').iso,tw/2) ] )	
	

	mdFile.new_header(level=2, title='Detector Operations')
	#mdFile.new_table(columns=)
		
	print(mpl.rcParams['text.usetex'], '3')

	f = GRBFollowup(name, skymap, start, stop, tw)
	f._allow_neg = args.allow_neg_ts
	f._dataset = 'GFUOnline_v001p02'
	f._fix_index = True
	f._float_index = False
	f._season_names = ['IC86, 2017', 'IC86, 2018', 'IC86, 2019']
	
	#f.initialize_injector() 
	#run_table = f.run_table(time_window=[1,2])
	
	run_table_list = ['Run', 'Start Time', 'Stop Time', 'Duration', 'Livetime']
	run_status_list = ['Run', 'Status', 'Light', 'Filter Mode', 'Run Mode', 'OK','GFU']
	length = 0

	run_url = 'https://live.icecube.wisc.edu/run_info/'
	run_query = {'user':'icecube', 'pass':'skua',
			'start':(Time(start_trigger[0]-0.33, format='mjd', precision=0)).iso,
			'stop': (Time(stop_trigger[0]+0.33, format = 'mjd',precision=0)).iso}
	now = Time(datetime.datetime.now()+datetime.timedelta(hours=6),scale='utc',precision=0)
	if use_urllib2:
		run_table = json.loads(urllib2.urlopen(
		urllib2.Request(run_url, urllib.urlencode(run_query)),
		timeout=500).read())
		print('this one')
	else:
		req_data = urllib.parse.urlencode(run_query).encode("utf-8")
		req = urllib.request.Request(run_url)
		with urllib.request.urlopen(req, data=req_data, timeout=500) as fi:
			run_table = json.load(fi)
	#add duration to run table

	ontime = {}
	ontime['type'] = 'database'
	ontime['stream'] = 'neutrino' ###NEED TO CHANGE FOR PREVIOUS MDJ TIMES####
	ontime['time_start'] = Time(run_table[0]['start'], format='iso', scale='utc', precision=0).iso
	ontime['time_stop'] = Time(run_table[-1]['stop'], format='iso', scale='utc', precision=0).iso

	query_events = icecube.realtime_tools.live.get_events(ontime['stream'], ontime['time_stop'], ontime['time_stop'])

	print(mpl.rcParams['text.usetex'], '4')

	'''
	widewindow = ontime_table(query_events)
	try:
		widewindow['t']=Time(list(widewindow['eventtime']),scale='utc',format='iso')
		for run in run_table:
			run['gfu_counts']=(widewindow['run_id'])
			run_status_list.append(run['gfu_counts'])
	except:
		print("Old years of data have different database keys")
		for run in run_table:
			run['gfu_counts']=0.
			run_status_list.append(run['gfu_counts'])
	'''

	for run in run_table:
		runstart = Time(run['start'],format='iso',scale='utc')
		if run['stop'] is None:
			runstop  = now
			run['stop']=runstop.iso
		else:
			runstop  = Time(run['stop'],scale='utc',precision=0)
		dt=int((runstop-runstart).sec)
		hours = int(dt)//3600
		minutes = int(dt)//60-hours*60
		seconds = dt%60

		livetime = (min(Time(stop_time,format='mjd'), runstop)-max(Time(start_time,format='mjd'), runstart)).sec
		run['livetime'] = livetime if livetime > 0 else 0
		if run['livetime'] > 86400. * 2:
			run['livetime'] = 0
			run['stop'] = run['start']
			dt=0
			hours =0
			minutes = 0
			seconds = 0

		run['duration_sec'] = dt
		run['duration'] = "{:01d}:{:02d}:{:02d}".format(hours,minutes,seconds)
		run['OK'] = ("OK" if run['status'] in ['SUCCESS']
			and run['lightmode'] in ['dark']
			and run['filter_mode'] in ['PhysicsFiltering']
			and run['run_mode'] in ['PhysicsTrig']
			else "NotOK"
			)
		
		#run['gfu_counts'] = (self.widewindow['run_id']==run['run_number']).sum()
		run_table_list.append(run['run_number'])
		run_table_list.append(run['start'])
		run_table_list.append(run['stop']) 
		run_table_list.append(run['duration'])
		run_table_list.append('{}s'.format(run['livetime']))
		
		run_status_list.append(run['run_number'])
		run_status_list.append(run['status'])
		run_status_list.append(run['lightmode'])
		run_status_list.append(run['filter_mode'])
		run_status_list.append(run['run_mode'])
		run_status_list.append(run['OK'])
		widewindow = ontime_table(query_events)
		try:
			widewindow['t']=Time(list(widewindow['eventtime']),scale='utc',format='iso')
			run['gfu_counts']=(widewindow['run_id'])
			run_status_list.append(run['gfu_counts'])
		except:
			print("Old years of data have different database keys")
			run['gfu_counts']=0.
			run_status_list.append(run['gfu_counts'])
		
		length = length + 1
	
	#print(run_table_list)
	mdFile.new_header(level=2, title='Run Times')
	mdFile.new_table(columns = 5, rows=length+1, text=run_table_list)
	mdFile.new_paragraph('Total Livetime = {}s'.format(run['livetime']))

	mdFile.new_header(level=2, title='Run Status')
	mdFile.new_table(columns = 7, rows=length+1, text=run_status_list)
	
	plot_path = '/data/user/rpmurphy/realtime/fast_response/FastResponseAnalysis/fast_response/scripts/2018_07_21_GRB180721A{}'.format(tw)
	'''make_rate_plots(start2stop, run_table, query_events, plot_path, args.name, tw, report_path,ontime['stream'])
	
	mdFile.new_header(level=2, title='Event Rates')
	mdFile.new_line(mdFile.new_inline_image(text='multiplicity', path='./IN_ICE_SIMPLE_MULTIPLICITY_plot_{}.png'.format(tw)))
	mdFile.new_line(mdFile.new_inline_image(text='badnessplot', path='./badness_plot_{}.png'.format(tw)))
	mdFile.new_line(mdFile.new_inline_image(text='muonfilter',path = './MuonFilter_13_plot_{}.png'.format(tw)))
	mdFile.new_line(mdFile.new_inline_image(text='Lfilter',path='./OnlineL2Filter_17_plot_{}.png'.format(tw)))
	mdFile.new_line(mdFile.new_inline_image(text='gfurate',path = './GFU_rate_plot.png_{}'.format(tw)))


	'''
	print('unblind TS')
	ts, ns = f.unblind_TS()
	print(mpl.rcParams['text.usetex'], '5')
	ts_list.append(ts)
	ns_list.append(ns)
	print('plot ontime')
	print(mpl.rcParams['text.usetex'], '6')
	f.plot_ontime()
	print('calc_pvalue')
	f.calc_pvalue()
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
	#y = f.upper_limit()
	#print(y)
	print('find coincident events')
	f.find_coincident_events()
	print('per event pvalue')
	f.per_event_pvalue()
	results = f.save_results()
	#f.generate_report()
	#print(plot_path)
	
	mdFile.new_header(level=1, title= 'Results')
	mdFile.new_line(mdFile.new_inline_image(text='All Sky On-time Events', path='./unblinded_skymap.png'))
	mdFile.new_line(mdFile.new_inline_image(text='Zoom Skymap', path='./unblinded_skymap_zoom.png'))
	mdFile.new_header(level=2, title = 'Coincident Events')
	#mdFile.new_line()
	### ADD COINCIDENT EVENT TABLE
	mdFile.new_header(level =2, title='Likelihood Analysis')
	mdFile.new_table(columns=4, rows=2, text=['ns', 'TS', 'p-value','sigma', ns, ts, p, sigma])
	mdFile.new_line(mdFile.new_inline_image(text='dN/dE', path = './central_90_dNdE_{}.png'.format(tw)))
	
	if os.path.exists(f'{plot_path}/TS_distribution_{tw}.png'):
		mdFile.new_line(mdFile.new_inline_image(text='TS distribution', path = './TS_distribution_{}.png'.format(tw)))
	else:
		mdFile.new_line('No background TS distrubution')

	if os.path.exists('./upper_limit_distribution.png'):
		mdFile.new_line(mdFile.new_inline_image(text='upper_limit_distribution.png'))
	else:
		mdFile.new_line('No upper limit calculation')

	if os.path.isdir(report_path):
        	print("Moving plots to {}".format(report_path))
       		subprocess.call(['cp {}/* {}'.format(plot_path, report_path)], shell=True)

	mdFile.new_line('These plots are meant as validations of the minimizer, and cannot be interpreted directly using Wilk’s assuming one degree of freedom because the low statistics invalidates the asymptotic condition for Wilk’s theorem.')
	
	np.delete(run_status_list, [0, len(run_status_list)-1])
	np.delete(run_table_list, [0, len(run_table_list)-1])
	#np.delete(ts_list)
	#np.delete(ns_list)

TW_best = np.argmin(p_vals)
p_post_best = p_vals[TW_best]
TS_unblind_best = ts_list[TW_best]



np.savez(report_path+f'{args.name}.npz',TS_unblind_list=ts_list, ns_unblind_list=ns_list, p_post=p_vals, TW_best=TW_best, p_post_best=p_post_best, TS_unblind_best=TS_unblind_best)

print('130427B_offline.npz'[:8])

mdFile.create_md_file()
text_file= open("GRB180721A_report.md","r")
text = text_file.read()
markdownFile = markdown.markdown(text, extensions=['tables', TableExtension(use_align_attribute=True)])
with open(report_path + "/{}_report.html".format(args.name, args.name), 'w') as output_file:
	output_file.write(markdownFile)


#markdown.markdownFromFile(input='markdownFile', output ='GRB180721A_report.html')

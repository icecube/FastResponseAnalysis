#!/usr/bin/env python

''' Script to automatically receive GCN notices for IceCube
    alert events and run followup accordingly

    Author: Alex Pizzuto, Jessie Thwaites, Alicia Mand
    Updated Date:   August 2026
'''

import logging
from gcn_kafka import Consumer
import os, subprocess, time, pwd, argparse
import healpy as hp
import numpy as np
import lxml.etree
from astropy.time import Time
from datetime import datetime
from dateutil.parser import parse
from glob import glob
from fast_response.slack_posters.slack import slackbot
import pandas as pd

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.warning("Connecting to GCN as Consumer")

with open('/home/jthwaites/private/tokens/kafka_token.txt') as f:
    client_id = f.readline().rstrip('\n')
    client_secret = f.readline().rstrip('\n')

domain = 'gcn.nasa.gov'
config = {'broker.address.family': 'v4', 
          'log_level': 0}
consumer = Consumer(client_id=client_id,
                    client_secret=client_secret,
                    domain='gcn.nasa.gov',
                    config=config,
                    #config={'max.poll.interval.ms':1800000},
                   )

#consumer.subscribe(['gcn.notices.icecube.gold_bronze_track_alerts'])
# stick with voevent for now for all 3
consumer.subscribe(['gcn.classic.voevent.ICECUBE_ASTROTRACK_BRONZE',
                    'gcn.classic.voevent.ICECUBE_ASTROTRACK_GOLD',
                    'gcn.classic.voevent.ICECUBE_CASCADE'])

def process_gcn(record): #payload,root
    analysis_path = os.environ.get('FAST_RESPONSE_SCRIPTS')
    if analysis_path is None:
        try:
            import fast_response
            analysis_path = os.path.dirname(fast_response.__file__) + '/scripts/'
        except Exception as e:
            logger.error('Error finding FRA package!!')
            post_error("Error finding FRA package")
            print('###########################################################################')
            print('CANNOT FIND ENVIRONMENT VARIABLE POINTING TO REALTIME FAST RESPONSE PACKAGE\n')
            print('You can either (1) install fast_response via pip or ')
            print('(2) put \'export FAST_RESPONSE_SCRIPTS=/path/to/fra/scripts\' in your bashrc')
            print('###########################################################################')
            raise Exception(e)

    # Read all of the VOEvent parameters from the "What" section.
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in record.iterfind('.//Param')}

    stream = params['Stream']
    eventtime = record.find('.//ISOTime').text
    if stream == '26':
        logger.warning("INCOMING ALERT: ",datetime.utcnow())
        logger.warning("Detected cascade type alert, running cascade followup. . . ")
        alert_type='cascade'
        event_name='IceCube-Cascade_{}{}{}'.format(eventtime[2:4],eventtime[5:7],eventtime[8:10])

        skymap = params['skymap_fits']
    else:
        alert_type='track'
        event_name='IceCube-{}{}{}'.format(eventtime[2:4],eventtime[5:7],eventtime[8:10]) 

        # IceCube sends 2: a notice and a revision, only want to run once
        if int(params['Rev']) !=0:
            return
        
        logger.warning("INCOMING ALERT: ",datetime.utcnow())
        logger.warning("Found track type alert, running track followup. . . ")

    event_id = params['event_id']
    run_id = params['run_id']
    event_mjd = Time(eventtime, format='isot').mjd

    # send message to slack with alert info
    bot = slackbot('fra-shifting')
    message =f'Listener found {alert_type} type alert, {event_name}. Waiting 1 day to run FRA'
    bot.post_short_msg(message)

    if alert_type == 'cascade':
        command = analysis_path + 'run_cascade_followup.py'
    else:
        command = analysis_path + 'run_track_followup.py'

    current_mjd = Time(datetime.utcnow(), scale='utc').mjd
    needed_delay = 1.
    current_delay = current_mjd - event_mjd
    while current_delay < needed_delay:
        logger.info("Need to wait another {:.1f} seconds before running".format(
            (needed_delay - current_delay)*86400.)
            )
        time.sleep((needed_delay - current_delay)*86400.)
        current_mjd = Time(datetime.utcnow(), scale='utc').mjd
        current_delay = current_mjd - event_mjd

    if alert_type == 'track':
        base_skymap_path = '/home/followup/output_plots/'
        skymap_f = glob(base_skymap_path \
            + f'run{int(run_id):08d}.evt{int(event_id):012d}.*probability.fits.gz')
        if len(skymap_f) == 0:
            logger.error("COULD NOT FIND THE SKYMAP FILE FOR V2 TRACK ALERT EVENT")
            return
        elif len(skymap_f) == 1:
            skymap = skymap_f[0]
        else:
            logger.error("TOO MANY OPTIONS FOR THE SKYMAP FILE FOR V2 TRACK ALERT EVENT")
            return
    
    #checking for events on the same day: looks for existing output files from previous runs
    count_dir=0
    for directory in os.listdir(os.environ.get('FAST_RESPONSE_OUTPUT')):
        if event_name in directory: count_dir+=1
    if count_dir==0: suffix='A'
    elif count_dir==2: suffix='B'
    elif count_dir==4: suffix='C'
    else: 
        logger.error("COULD NOT DETERMINE EVENT SUFFIX")
        logger.error("check for other events on the same day and re-run with args:")
        logger.error('--skymap={} --time={} --alert_id={}'.format(skymap, str(event_mjd), run_id+':'+event_id)) 
        return
    
    logger.info('\nRunning {} --skymap={} --time={} --alert_id={} --suffix={}'.format(
        command, skymap, str(event_mjd), run_id+':'+event_id, suffix))

    print(params)

    subprocess.call([command, '--skymap={}'.format(skymap), 
        '--time={}'.format(str(event_mjd)), 
        '--alert_id={}'.format(run_id+':'+event_id),
        '--suffix={}'.format(suffix)]
        )

    event_name=event_name+suffix
    doc = False
    if args.document:
        try:
            dir_1000 = glob(os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'),
                                          '*{}_1.0e+03_s').format(event_name))
            subprocess.call([analysis_path+'document.py', '--path', dir_1000[0]])
            dir_2d = glob(os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'),
                                          '*{}_1.7e+05_s').format(event_name))
            subprocess.call([analysis_path+'document.py', '--path', dir_2d[0]])
            doc=True
        except:
            post_error("Failed to run document command")
            logger.warning('Failed to document to private webpage')

    try: 
        shifters = pd.read_csv(os.path.join(analysis_path,'../slack_posters/fra_shifters.csv'), 
                               parse_dates=[0,1])
        on_shift=''
        for i in shifters.index:
            if shifters['start'][i]<datetime.utcnow()<shifters['stop'][i]:
                on_shift+='<@{}> '.format(shifters['slack_id'][i])
        link = 'https://user-web.icecube.wisc.edu/~jthwaites/FastResponse/webpage/output/'
        start_1000 = Time(event_mjd -500./86400., format='mjd').iso
        wp_link_1000 = '{}{}_{}_1.0e+03_s.html'.format(link, start_1000[:10].replace('-','_'),event_name)
        
        start_2d = Time(event_mjd-1., format='mjd').iso
        wp_link_2d   = '{}{}_{}_1.7e+05_s.html'.format(link, start_2d[:10].replace('-','_'), event_name)
        done_message = f'Done running FRA for {alert_type} alert, {event_name}.\n '+ on_shift +'on shift'

        if doc:
            done_message = done_message + "\n - Results for 1000s: <{}|link>\n - Results for 2d: <{}|link>".format(
                              wp_link_1000, wp_link_2d)

        bot.post_short_msg(done_message)
    except Exception as e:
        post_error("Failed to push results to private webpage")
        logger.warning('Failed to push to private webpage')
        logger.warning(e)

def post_error(errMsg=None): 
    analysis_path = os.environ.get('FAST_RESPONSE_SCRIPTS')
    bot = slackbot('fra-shifting')
    if errMsg != None: 
        message = errMsg
    else: 
        message = "ERROR in FRA, please check internal listener!"
    try: 
        shifters = pd.read_csv(os.path.join(analysis_path, '../slack_posters/fra_shifters.csv'), parse_dates=[0,1])
        on_shift = ''
        for i in shifters.index:
            if shifters['start'][i] < datetime.utcnow() < shifters['stop'][i]: 
                on_shift+='<@{}> '.format(shifters['slack_id'][i])
        error_message = f"{message} {on_shift} on shift."
        bot.post_short_msg(error_message)
    except Exception as e: 
        logger.warning("Failed to post error message")
        logger.warning(e)
    return 

if __name__ == '__main__':

    username = pwd.getpwuid(os.getuid())[0]

    #default for if to document or not: only way to check reports on realtime 
    if username == 'realtime': 
        document = True
    else: 
        document = False

    parser = argparse.ArgumentParser(description='Fast Response Analysis')
    parser.add_argument('--run_live', action='store_true', default=False,
                        help='Run on live GCNs')
    parser.add_argument('--test_cascade', default=False, action='store_true',
                        help='When testing, raise to run a cascade, else track')
    parser.add_argument('--test_error', default=False, action='store_true', 
                        help='When testing, raise to post an error message')
    parser.add_argument('--document', action='store_true', default=document,
                        help='flag to raise to push results to internal webpage')
    args = parser.parse_args()

    if args.run_live:
        logger.warning("Listening for IC Alert GCNs . . . ")
        while True:
            for message in consumer.consume(timeout=1):
                if message.error():
                    logger.warning(message.error())
                    continue
                # value = message.value().decode('utf-8')
                # value = value.replace("<?xml version='1.0' encoding='UTF-8'?>","") #lxml doesn't like this line
                value = message.value()
                logger.warning('Found GCN on topic {}'.format(message.topic()))
                notice = lxml.etree.fromstring(value)
                try: 
                    process_gcn(notice)
                except Exception as e:
                    post_error("Could not process GCN")
                    logger.warning("Could not process GCN: ", e) 

    else:
        try:
            import fast_response
            sample_skymap_path=os.path.dirname(fast_response.__file__) +'/sample_skymaps/'
        except Exception as e:
            post_error("Cannot find path to sample skymaps")
            logger.error('Cannot find path to sample skymaps')
            raise Exception(e)
        if args.test_error: 
            post_error()

        if not args.test_cascade and not args.test_error:
            logger.info("Running on sample track . . . ")
            payload = open(sample_skymap_path \
                + 'sample_astrotrack_alert_2021.xml', 'rb').read()
            root = lxml.etree.fromstring(payload)
            process_gcn(root)
        elif args.test_cascade and not args.test_error:
            logger.info("Running on sample cascade . . . ")
            payload = open(sample_skymap_path \
                + 'sample_cascade.txt', 'rb').read()
            root = lxml.etree.fromstring(payload)
            process_gcn(root)

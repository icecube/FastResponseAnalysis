#!/usr/bin/env python

#rm ln 103/4!!!
import logging
from datetime import datetime
import socket
import requests
import healpy as hp
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import io, time, os
import urllib.request, urllib.error, urllib.parse
import argparse
import json, pickle
from gcn_kafka import Consumer
from icecube import realtime_tools
import numpy as np
import lxml.etree
from astropy.time import Time
import dateutil.parser
from datetime import datetime

with open('/cvmfs/icecube.opensciencegrid.org/users/jthwaites/tokens/kafka_token.txt') as f:
    client_id = f.readline().rstrip('\n')
    client_secret = f.readline().rstrip('\n')

consumer = Consumer(client_id=client_id,
                    client_secret=client_secret)

# Subscribe to topics to receive alerts
consumer.subscribe(['gcn.classic.voevent.LVC_PRELIMINARY',
                    'gcn.classic.voevent.LVC_INITIAL',
                    'gcn.classic.voevent.LVC_UPDATE'])

def SendAlert(results=None):
        from gcn_kafka import Producer

        if results is None:
            logger.fatal('Found no alert to send')
        
        with open('/cvmfs/icecube.opensciencegrid.org/users/jthwaites/tokens/real_icecube_kafka_prod.txt') as f:
            prod_id = f.readline().rstrip('\n')
            prod_secret = f.readline().rstrip('\n')

        producer = Producer(client_id=prod_id,
                            client_secret=prod_secret,
                            domain='gcn.nasa.gov')
        
        if 'MS' in results['ref_id']:
            return
            topic = 'gcn.notices.icecube.test.lvk_nu_track_search'
        else:
            topic = 'gcn.notices.icecube.lvk_nu_track_search'
        
        logger.info('sending to {} on gcn.nasa.gov'.format(topic))
        sendme = json.dumps(results)
        producer.produce(topic, sendme.encode())
        ret = producer.flush()
        return ret

def SendTestAlert(results=None):
        from gcn_kafka import Producer

        if results is None:
            logger.fatal('Found no alert to send')
        
        with open('/cvmfs/icecube.opensciencegrid.org/users/jthwaites/tokens/test_icecube_kafka_prod.txt') as f:
            prod_id = f.readline().rstrip('\n')
            prod_secret = f.readline().rstrip('\n')

        producer = Producer(client_id=prod_id,
                            client_secret=prod_secret,
                            domain='test.gcn.nasa.gov')
        
        topic = 'gcn.notices.icecube.test.lvk_nu_track_search'
        logger.info('sending to {} on test.gcn.nasa.gov'.format(topic))
        
        sendme = json.dumps(results)
        producer.produce(topic, sendme.encode())
        ret = producer.flush()
        return ret

def format_ontime_events_uml(events, event_mjd):
    ontime_events={}
    for event in events:
        ontime_events[event['event']]={
            'event_dt' : round((event['time']-event_mjd)*86400., 2),
            'localization':{
                'ra' : round(np.rad2deg(event['ra']), 2),
                'dec' : round(np.rad2deg(event['dec']), 2),
                "uncertainty_shape": "circle",
                'ra_uncertainty': round(np.rad2deg(event['sigma']*2.145966),2),
                "containment_probability": 0.9,
                "systematic_included": False
            },
            'event_pval_generic' : round(event['pvalue'],4)
        }
    return ontime_events

def format_ontime_events_llama(events):
    ontime_events={}
    for event in events:
        ontime_events[event['i3event']] = {
            'event_dt' : round(event['dt'],2),
            'localization':{
                'ra' : round(event['ra'], 2),
                'dec' : round(event['dec'], 2),
                "uncertainty_shape": "circle",
                'ra_uncertainty': round(np.rad2deg(np.deg2rad(event['sigma'])*2.145966),3),
                "containment_probability": 0.9,
                "systematic_included": False
            },
            'event_pval_bayesian': round(event['p_value'],4)
        }
    return ontime_events

def format_ontime_events_llama_old(events,event_mjd):
    ontime_events={}
    for event in events:
        ontime_events[event['event']] = {
            'event_dt' : round((event['mjd']-event_mjd)*86400.,2),
            'localization':{
                'ra' : round(np.rad2deg(event['ra']), 2),
                'dec' : round(np.rad2deg(event['dec']), 2),
                "uncertainty_shape": "circle",
                'ra_uncertainty': round(np.rad2deg(event['sigma']*2.145966),3),
                "containment_probability": 0.9,
                "systematic_included": False
            },
            'event_pval_bayesian': 'null'
        }
    return ontime_events

def combine_events(uml_ontime, llama_ontime):
    #in the case both return coincident events, want to combine their results
    event_ids_all = np.unique(list(uml_ontime.keys())+list(llama_ontime.keys()))
    coinc_events = []

    for id in event_ids_all:
        #case: event in both results
        if (id in uml_ontime.keys()) and (id in llama_ontime.keys()):
            if (uml_ontime[id]['event_pval_generic'] < 0.1) and (llama_ontime[id]['event_pval_bayesian']<0.1):
                uml_ontime[id]['event_pval_bayesian'] = llama_ontime[id]['event_pval_bayesian']
                coinc_events.append(uml_ontime[id])
            elif uml_ontime[id]['event_pval_generic'] < 0.1:
                uml_ontime[id]['event_pval_bayesian'] = 'null'
                coinc_events.append(uml_ontime[id])
            elif llama_ontime[id]['event_pval_bayesian']<0.1:
                uml_ontime[id]['event_pval_bayesian'] = llama_ontime[id]['event_pval_bayesian']
                uml_ontime[id]['event_pval_generic'] ='null'
                coinc_events.append(uml_ontime[id])
        #case: only in uml
        elif id in uml_ontime.keys():
            if uml_ontime[id]['event_pval_generic'] < 0.1:
                uml_ontime[id]['event_pval_bayesian'] = 'null'
                coinc_events.append(uml_ontime[id])
        #case: only in llama
        elif id in llama_ontime.keys():
            if llama_ontime[id]['event_pval_bayesian']<0.1:
                llama_ontime[id]['event_pval_generic']='null'
                coinc_events.append(llama_ontime[id])

    return coinc_events

def parse_notice(record, wait_for_llama=False, heartbeat=False):
    logger = logging.getLogger()

    if record.attrib['role']!='observation':
        logger.warning('found test event')
        fra_results_location = '/data/user/jthwaites/o4-mocks/'
        if not heartbeat:
            return
    else:
        logger.warning('ALERT FOUND')
        fra_results_location = os.environ.get('FAST_RESPONSE_OUTPUT')
    
    ### GENERAL PARAMETERS ###
    #read event information
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in record.iterfind('.//Param')}

    # ignore subthreshold, only run on significant events
    subthreshold=False
    if 'Significant' in params.keys():
        if int(params['Significant'])==0: 
            subthreshold=True
    if params['Group'] == 'Burst':
        wait_for_llama = False
        m = 'Significant' if not subthreshold else 'Subthreshold'
        logger.warning('{} burst alert found. '.format(m))

    if not wait_for_llama and subthreshold:
        #if llama isn't running, don't run on subthreshold
        logger.warning('Not waiting for LLAMA, subthreshold alert. Returning...')
        return

    logger.info("alert found, processing GCN")

    collected_results = {}
    collected_results["$schema"]= "https://gcn.nasa.gov/schema/gcn/notices/icecube/LvkNuTrackSearch.schema.json"
    collected_results["type"]= "IceCube LVK Alert Nu Track Search"

    eventtime = record.find('.//ISOTime').text
    event_mjd = Time(eventtime, format='isot').mjd
    name = record.attrib['ivorn'].split('#')[1] 

    if record.attrib['role'] != 'observation':
        name=name+'_test'

    collected_results['ref_id'] = name

    #collected_results['ref_id'] = name.split('-')[0]
    #collected_results['id'] = [name.split('-')[2],'Sequence:{}'.format(name.split('-')[1])]
    
    ### WAIT ###
    # wait until the 500s has elapsed for data
    current_mjd = Time(datetime.utcnow(), scale='utc').mjd
    needed_delay = 500./84600. 
    current_delay = current_mjd - event_mjd

    while current_delay < needed_delay:
        logger.info("Wait {:.1f} seconds for data".format(
            (needed_delay - current_delay)*86400.)
            )
        time.sleep((needed_delay - current_delay)*86400.)
        current_mjd = Time(datetime.utcnow(), scale='utc').mjd
        current_delay = current_mjd - event_mjd

    # start checking for results
    results_done = False
    start_check_mjd = Time(datetime.utcnow(), scale='utc').mjd
    max_delay = max_wait/60./24. #mins to days

    while results_done == False:
        logger.info("Waiting for results (max {:.0f} mins)".format(
            max_wait)
            )
        start_date = Time(dateutil.parser.parse(eventtime)).datetime
        start_str = f'{start_date.year:02d}_{start_date.month:02d}_{start_date.day:02d}'

        uml_results_path = os.path.join(fra_results_location, start_str + '_' + name.replace(' ', '_') \
                                   + '/' + start_str + '_' + name.replace(' ', '_')+'_results.pickle')
        uml_results_finished = os.path.exists(uml_results_path)

        llama_name = '{}.significance_subthreshold_lvc-i3.json'.format(record.attrib['ivorn'].split('#')[1])
        llama_results_path = os.path.join(llama_results_location,llama_name)
        if wait_for_llama:
            llama_results_finished = os.path.exists(llama_results_path)
        else: 
            llama_results_finished = False

        if subthreshold and llama_results_finished:
            #no UML results in subthreshold case, only LLAMA
            results_done=True
            logger.info('found results for LLAMA for subthreshold event, writing notice')

        if uml_results_finished and llama_results_finished:
            results_done=True
            logger.info('found results for both, writing notice')
        elif uml_results_finished and not wait_for_llama:
            results_done = True
            logger.info('found results for UML, writing notice')
        else:
            current_mjd = Time(datetime.utcnow(), scale='utc').mjd
            if current_mjd - start_check_mjd < max_delay:
                time.sleep(10.)
            else:
                if (uml_results_finished or llama_results_finished):
                    if uml_results_finished:
                        logger.warning('LLAMA results not finished in {:.0f} mins. Sending UML only'.format(max_wait))
                        results_done=True
                    if llama_results_finished:
                        logger.warning('UML results not finished in {:.0f} mins. Sending LLAMA only'.format(max_wait))
                        results_done=True
                else:
                    logger.warning('Both analyses not finished after {:.0f} min wait.'.format(max_wait))
                    logger.warning('Not sending GCN.')
                    return
                
    collected_results['alert_datetime'] = '{}Z'.format(Time(datetime.utcnow(), scale='utc').isot)
    collected_results['trigger_time'] = eventtime if 'Z' in eventtime else '{}Z'.format(eventtime) 
    collected_results['observation_start'] = '{}Z'.format(Time(event_mjd - 500./86400., format = 'mjd').isot)
    collected_results['observation_stop'] = '{}Z'.format(Time(event_mjd + 500./86400., format = 'mjd').isot)
    collected_results['observation_livetime'] = 1000

    ### COLLECT RESULTS ###
    additional_website_params = {}
    if uml_results_finished:
        with open(uml_results_path, 'rb') as f:
            uml_results = pickle.load(f)

        for key in ['skymap_path', 'analysisid']:
            if key in uml_results.keys():
                additional_website_params[key] = uml_results[key]

    if llama_results_finished:
        with open(llama_results_path, 'r') as f:
            llama_results = json.load(f)
        if (record.attrib['role'] == 'observation') and (llama_results['inputs']['neutrino_info'][0]['type']=='blinded'):
            logger.warning('LLAMA results blinded for real event! Skipping LLAMA')
            llama_results_finished = False
            if subthreshold: 
                logger.warning('LLAMA skipped on subthreshold event, returning')
                return

    if uml_results_finished and llama_results_finished:
        collected_results['pval_generic'] = round(uml_results['p'],4)
        collected_results['pval_bayesian'] = round(llama_results['p_value'],4)
        if (collected_results['pval_generic']<0.1) or (collected_results['pval_bayesian']<0.1):
            uml_ontime = format_ontime_events_uml(uml_results['coincident_events'], event_mjd)
            llama_ontime = format_ontime_events_llama(llama_results['single_neutrino'])
            coinc_events = combine_events(uml_ontime, llama_ontime)
            collected_results['n_events_coincident'] = len(coinc_events)
            collected_results['coincident_events'] = coinc_events

            collected_results['most_likely_direction'] = {
                'ra': round(np.rad2deg(uml_results['fit_ra']), 2),
                'dec' : round(np.rad2deg(uml_results['fit_dec']), 2),
            }
        else: 
            collected_results['n_events_coincident'] = 0
        
        collected_results["neutrino_flux_sensitivity_range"] = {
                'flux_sensitivity' : [
                    round(uml_results['sens_range'][0],4),
                    round(uml_results['sens_range'][1],4)
                ],
                'sensitive_energy_range' : [
                    round(float('{:.2e}'.format(uml_results['energy_range'][0]))),
                    round(float('{:.2e}'.format(uml_results['energy_range'][1])))
                ],
        }
    
    elif uml_results_finished:
        collected_results['pval_generic'] = round(uml_results['p'],4)
        collected_results['pval_bayesian'] = 'null'
        if collected_results['pval_generic'] <0.1:
            uml_ontime = format_ontime_events_uml(uml_results['coincident_events'], event_mjd)
            coinc_events=[]
            for eventid in uml_ontime.keys():
                if (uml_ontime[eventid]['event_pval_generic'] < 0.1):
                    uml_ontime[eventid]['event_pval_bayesian'] = 'null'
                    coinc_events.append(uml_ontime[eventid])
            collected_results['n_events_coincident'] = len(coinc_events)
            collected_results['coincident_events'] = coinc_events

            collected_results['most_likely_direction'] = {
                'ra': round(np.rad2deg(uml_results['fit_ra']), 2),
                'dec' : round(np.rad2deg(uml_results['fit_dec']), 2)
            }

        else: 
            collected_results['n_events_coincident'] = 0
        
        collected_results["neutrino_flux_sensitivity_range"] = {
                'flux_sensitivity' : [
                    round(uml_results['sens_range'][0],4),
                    round(uml_results['sens_range'][1],4)
                ],
                'sensitive_energy_range' :[
                    round(float('{:.2e}'.format(uml_results['energy_range'][0]))),
                    round(float('{:.2e}'.format(uml_results['energy_range'][1])))
                ],
        }

    elif llama_results_finished:
        collected_results['pval_generic'] = 'null'
        collected_results['pval_bayesian'] = round(llama_results['p_value'],4)
        if collected_results['pval_bayesian']<0.1:
            llama_ontime = format_ontime_events_llama(llama_results['single_neutrino'])
            coinc_events=[]
            for eventid in llama_ontime.keys():
                if (llama_ontime[eventid]['event_pval_bayesian'] < 0.1):
                    llama_ontime[eventid]['event_pval_generic'] = 'null'
                    coinc_events.append(llama_ontime[eventid])
            collected_results['n_events_coincident'] = len(coinc_events)
            collected_results['coincident_events'] = coinc_events
        else: 
            collected_results['n_events_coincident'] = 0

        collected_results["neutrino_flux_sensitivity_range"] = {
                'flux_sensitivity' : [
                    round(llama_results['neutrino_flux_sensitivity_range']['flux_sensitivity'][0],4),
                    round(llama_results['neutrino_flux_sensitivity_range']['flux_sensitivity'][1],4)
                ],
                'sensitive_energy_range' :[
                    round(float('{:.2e}'.format(llama_results['neutrino_flux_sensitivity_range']
                                                ['sensitive_energy_range'][0]))),
                    round(float('{:.2e}'.format(llama_results['neutrino_flux_sensitivity_range']
                                                ['sensitive_energy_range'][1])))
                ],
        }

    if (collected_results['n_events_coincident'] == 0) and ('coincident_events' in collected_results.keys()):
        c = collected_results.pop('coincident_events')
    if ('most_likely_direction' in collected_results.keys()):
        try: 
            if (collected_results['pval_generic']>0.1):
                c = collected_results.pop('most_likely_direction')
        except Exception as e:
            print(e)
        
    ### SAVE RESULTS ###
    if record.attrib['role']=='observation' and not heartbeat:
        with open(os.path.join(save_location, f'{name}_collected_results.json'),'w') as f:
            json.dump(collected_results, f, indent = 6)

        logger.info('sending notice')
        status = SendAlert(results = collected_results)
        st = 'sent' if status == 0 else 'error!'
        logger.info('status: {}'.format(status))
        logger.info('{}'.format(st))
        
    else:
        with open(os.path.join(save_location, f'mocks/{name}_collected_results.json'),'w') as f:
            json.dump(collected_results, f, indent = 6)
        #status = SendTestAlert(results = collected_results)
        #logger.info('status: {}'.format(status))

    # Adding a few extra keys needed to create public webpage
    for key in additional_website_params.keys():
        collected_results[key] = additional_website_params[key]
    
    #web_utils.updateGW_public(collected_results)

parser = argparse.ArgumentParser(description='Combine GW-Nu results')
parser.add_argument('--run_live', action='store_true', default=False,
                    help='Run on live GCNs')
parser.add_argument('--test_path', type=str, default=None,
                    help='path to test xml file')
parser.add_argument('--max_wait', type=float, default=60.,
                    help='Maximum minutes to wait for LLAMA/UML results before timeout (default=60)')
parser.add_argument('--wait_for_llama', action='store_true', default=False,
                    help='bool to decide to send llama results with uml, default false while not unblinded')
parser.add_argument('--heartbeat', action='store_true', default=False,
                    help='bool to save jsons for mocks (saves to save_dir+mocks/)')
parser.add_argument('--save_dir', type=str, default='/home/followup/lvk_followup_output/',
                    help='Directory to save output json (default=/home/followup/lvk_followup_output/)')
args = parser.parse_args()

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.warning("combine_results starting, connecting to GCN")

#fra_results_location = os.environ.get('FAST_RESPONSE_OUTPUT')#'/data/user/jthwaites/o4-mocks/'
llama_results_location = '/home/followup/lvk_dropbox/'
#save_location = '/home/followup/lvk_followup_output/' #where to save final json
save_location = args.save_dir

max_wait = args.max_wait
#wait_for_llama = args.wait_for_llama

if args.run_live:
    logger.info('running on live GCNs')
    while True:
        for message in consumer.consume(timeout=1):
            value = message.value()
            if '<' not in value.decode('utf-8')[0]:
                #sometimes, we'll get error messages - these make the code fail. skip them
                logger.warning(value.decode('utf-8'))
                continue
            notice = lxml.etree.fromstring(value.decode('utf-8').encode('ascii'))
            parse_notice(notice, wait_for_llama=args.wait_for_llama, heartbeat = args.heartbeat)
            logger.info('Done.')

else:
    if args.test_path is None:
        import glob
        paths=glob.glob('/home/jthwaites/FastResponse/*/*xml')
        path = paths[1]
        #path = '/home/jthwaites/FastResponse/S230522n-preliminary.json,1'
    else: 
        path = args.test_path
    
    logger.info('running on {}'.format(path))
    
    with open(path, 'r') as f:
        payload = f.read()
    try:
        record = lxml.etree.fromstring(payload)
    except Exception as e:
        print(e)
        exit()
    
    parse_notice(record, wait_for_llama=args.wait_for_llama)
logger.info("done")
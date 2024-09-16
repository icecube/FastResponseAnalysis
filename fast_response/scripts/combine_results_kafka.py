#!/usr/bin/env python

import logging
from datetime import datetime
import socket
import requests
import healpy as hp
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import io, time, os, glob, subprocess
import urllib.request, urllib.error, urllib.parse
import argparse
import json, pickle
from gcn_kafka import Consumer
#from icecube import realtime_tools
import numpy as np
import lxml.etree
from astropy.time import Time
import dateutil.parser
from datetime import datetime
from fast_response.web_utils import updateGW_public

with open('/home/jthwaites/private/tokens/kafka_token.txt') as f:
    client_id = f.readline().rstrip('\n')
    client_secret = f.readline().rstrip('\n')

consumer = Consumer(client_id=client_id,
                    client_secret=client_secret,
                    config={'max.poll.interval.ms':1800000})

# Subscribe to topics to receive alerts
consumer.subscribe(['gcn.classic.voevent.LVC_PRELIMINARY',
                    'gcn.classic.voevent.LVC_INITIAL',
                    'gcn.classic.voevent.LVC_UPDATE'])
#consumer.subscribe(['igwn.gwalert'])

def SendAlert(results=None):
        from gcn_kafka import Producer

        if results is None:
            logger.fatal('Found no alert to send')
        
        with open('/home/jthwaites/private/tokens/real_icecube_kafka_prod.txt') as f:
            prod_id = f.readline().rstrip('\n')
            prod_secret = f.readline().rstrip('\n')

        producer = Producer(#config={'bootstrap.servers': 'kafka3.gcn.nasa.gov'},
                            client_id=prod_id,
                            client_secret=prod_secret,
                            domain='gcn.nasa.gov')
        
        if 'MS' in results['ref_ID']:
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
        
        with open('/home/jthwaites/private/tokens/test_icecube_kafka_prod.txt') as f:
            prod_id = f.readline().rstrip('\n')
            prod_secret = f.readline().rstrip('\n')

        producer = Producer(#config={'bootstrap.servers': 'kafka3.gcn.nasa.gov'},
                            client_id=prod_id,
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
                'ra_dec_error': round(np.rad2deg(event['sigma']*2.145966),2),
                "containment_probability": 0.9,
                "systematic_included": False
            },
            #'event_pval_generic' : round(event['pvalue'],4)
        }
        if event['pvalue'] < 0.0001:
            ontime_events[event['event']]['event_pval_generic'] = float('{:.1e}'.format(event['pvalue']))
        else:
            ontime_events[event['event']]['event_pval_generic'] = round(event['pvalue'],4)
    return ontime_events

def format_ontime_events_llama(events):
    ontime_events={}
    for event in events:
        ontime_events[event['i3event']] = {
            'event_dt' : round(event['dt'],2),
            'localization':{
                'ra' : round(event['ra'], 2),
                'dec' : round(event['dec'], 2),
                'ra_dec_error': round(np.rad2deg(np.deg2rad(event['sigma'])*2.145966),2),
                "containment_probability": 0.9,
                "systematic_included": False
            },
            #'event_pval_bayesian': round(event['p_value'],4)
        }
        if event['p_value'] < 0.0001:
            ontime_events[event['i3event']]['event_pval_bayesian'] = float('{:.1e}'.format(event['p_value']))
        else:
            ontime_events[event['i3event']]['event_pval_bayesian'] = round(event['p_value'],4)
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
                uml_ontime[id]['event_pval_bayesian'] = None
                coinc_events.append(uml_ontime[id])
            elif llama_ontime[id]['event_pval_bayesian']<0.1:
                uml_ontime[id]['event_pval_bayesian'] = llama_ontime[id]['event_pval_bayesian']
                uml_ontime[id]['event_pval_generic'] = None
                coinc_events.append(uml_ontime[id])
        #case: only in uml
        elif id in uml_ontime.keys():
            if uml_ontime[id]['event_pval_generic'] < 0.1:
                uml_ontime[id]['event_pval_bayesian'] = None
                coinc_events.append(uml_ontime[id])
        #case: only in llama
        elif id in llama_ontime.keys():
            if llama_ontime[id]['event_pval_bayesian']<0.1:
                llama_ontime[id]['event_pval_generic']= None
                coinc_events.append(llama_ontime[id])

    return coinc_events

def parse_notice(record, wait_for_llama=False, heartbeat=False):
    logger = logging.getLogger()

    if record.attrib['role']!='observation':
        fra_results_location = '/data/user/jthwaites/o4-mocks/'
        if not heartbeat:
            logger.warning('found test event - not in mock mode. returning')
            return
        else:
            logger.warning('found test event')
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
            logger.warning('low-significance alert found. ')
    if params['Group'] == 'Burst' or params["Pipeline"] =='CWB':
        wait_for_llama = False
        m = 'Significant' if not subthreshold else 'Subthreshold'
        logger.warning('{} burst or CWB alert found. '.format(m))
    if len(params['Instruments'].split(','))==1:
        #wait_for_llama = False
        logger.warning('One detector event found. ')

    if not wait_for_llama and subthreshold:
        #if llama isn't running, don't run on subthreshold
        logger.warning('Not waiting for LLAMA. Returning...')
        return

    collected_results = {}
    collected_results["$schema"]= "https://gcn.nasa.gov/schema/v4.1.0/gcn/notices/icecube/lvk_nu_track_search.schema.json"
    collected_results["type"]= "IceCube LVK Alert Nu Track Search"

    eventtime = record.find('.//ISOTime').text
    event_mjd = Time(eventtime, format='isot').mjd
    name = record.attrib['ivorn'].split('#')[1] 

    if record.attrib['role'] != 'observation':
        name=name+'_test'

    logger.info("{} alert found, processing GCN".format(name))
    collected_results['ref_ID'] = name.split('-')[0]
    collected_results['reference']= {"gcn.notices.LVK.alert": name}
    
    ### WAIT ###
    # wait until the 500s has elapsed for data
    current_mjd = Time(datetime.utcnow(), scale='utc').mjd
    needed_delay = 500./84600. 
    #needed_delay = 720./86400. #12 mins while lvk_dropbox is offline
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

    logger.info("Looking for results (max {:.0f} mins)".format(
            max_wait)
            )
    
    while results_done == False:
        start_date = Time(dateutil.parser.parse(eventtime)).datetime
        start_str = f'{start_date.year:02d}_{start_date.month:02d}_{start_date.day:02d}'

        uml_results_path = os.path.join(fra_results_location, start_str + '_' + name.replace(' ', '_') \
                                   + '/' + start_str + '_' + name.replace(' ', '_')+'_results.pickle')
        uml_results_finished = os.path.exists(uml_results_path)

        if len(params['Instruments'].split(','))==1:
            llama_name = '{}.significance_opa_lvc-i3.json'.format(record.attrib['ivorn'].split('#')[1])
        else:
            llama_name = '{}.significance_subthreshold_lvc-i3.json'.format(record.attrib['ivorn'].split('#')[1])
        llama_results_path = os.path.join(llama_results_location,llama_name)
        if wait_for_llama:
            llama_results_finished = os.path.exists(llama_results_path)
        else: 
            llama_results_finished = False
        #llama_name = '{}*-significance_subthreshold_lvc-i3.json'.format(record.attrib['ivorn'].split('#')[1])
        #llama_results_path = os.path.join(llama_results_location,llama_name)
        #if wait_for_llama:
        #    llama_results_glob = sorted(glob.glob(llama_results_path))
        #    if len(llama_results_glob)>0:
        #        llama_results_path = llama_results_glob[-1]
        #        llama_results_finished=True
        #    else:
        #        llama_results_finished=False
        #else: 
        #    llama_results_finished = False

        if subthreshold and llama_results_finished:
            #no UML results in subthreshold case, only LLAMA
            results_done = True
            uml_results_finished = False
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
                        if record.attrib['role']=='observation' and not heartbeat:
                            try: 
                                subprocess.call(['/home/jthwaites/private/make_call.py', 
                                                 '--troubleshoot_gcn=True', '--missing_llama=True'])
                            except:
                                logger.warning('Failed to send alert to shifters: Issue finding LLAMA results. ')
                    if llama_results_finished:
                        logger.warning('UML results not finished in {:.0f} mins. Sending LLAMA only'.format(max_wait))
                        results_done=True
                        if record.attrib['role']=='observation' and not heartbeat:
                            try: 
                                subprocess.call(['/home/jthwaites/private/make_call.py', 
                                                 '--troubleshoot_gcn=True', '--missing_uml=True'])
                            except:
                                logger.warning('Failed to send alert to shifters: Issue finding UML results. ')
                else:
                    logger.warning('Both analyses not finished after {:.0f} min wait.'.format(max_wait))
                    logger.warning('Not sending GCN.')

                    if record.attrib['role']=='observation' and not heartbeat:
                        err_msg = ['/home/jthwaites/private/make_call.py', '--troubleshoot_gcn=True', 
                                   '--missing_llama=True']
                        if not subthreshold: err_msg.append('--missing_uml=True')
                        
                        try: 
                            subprocess.call(err_msg)
                        except:
                            logger.warning('Failed to send alert to shifters: Issue finding both results. ')
                    return
                
    collected_results['alert_datetime'] = '{}Z'.format(Time(datetime.utcnow(), scale='utc').isot)
    collected_results['trigger_time'] = eventtime if 'Z' in eventtime else '{}Z'.format(eventtime) 
    collected_results['observation_start'] = '{}Z'.format(Time(event_mjd - 500./86400., format = 'mjd').isot)
    collected_results['observation_stop'] = '{}Z'.format(Time(event_mjd + 500./86400., format = 'mjd').isot)
    collected_results['observation_livetime'] = 1000

    ### COLLECT RESULTS ###
    send_notif = False

    if uml_results_finished:
        with open(uml_results_path, 'rb') as f:
            uml_results = pickle.load(f)

    if llama_results_finished:
        with open(llama_results_path, 'r') as f:
            llama_results = json.load(f)
        try:
            if (record.attrib['role'] == 'observation') and (llama_results['inputs']['neutrino_info'][0]['type']=='blinded'):
                logger.warning('LLAMA results blinded for real event! Skipping LLAMA')
                llama_results_finished = False
                if subthreshold: 
                    logger.warning('LLAMA skipped on subthreshold event, returning')
                    return
        except:
            logger.warning('NO neutrinos in 1000s in LLAMA result! Still sending . . .')

    if uml_results_finished and llama_results_finished:
        if uml_results['p']< 0.0001:
            collected_results['pval_generic'] = float('{:.1e}'.format(uml_results['p']))
        else:
            collected_results['pval_generic'] = round(uml_results['p'],4)
        
        if llama_results['p_value'] < 0.0001:
            collected_results['pval_bayesian'] = float('{:.1e}'.format(llama_results['p_value']))
        else:
            collected_results['pval_bayesian'] = round(llama_results['p_value'],4)

        if (collected_results['pval_generic']<0.01) or (collected_results['pval_bayesian']<0.01):
            send_notif=True

        if (collected_results['pval_generic']<0.1) or (collected_results['pval_bayesian']<0.1):
            uml_ontime = format_ontime_events_uml(uml_results['coincident_events'], event_mjd)
            llama_ontime = format_ontime_events_llama(llama_results['single_neutrino'])
            coinc_events = combine_events(uml_ontime, llama_ontime)
            collected_results['n_events_coincident'] = len(coinc_events)
            collected_results['coincident_events'] = coinc_events

            collected_results['most_probable_direction'] = {
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
        if uml_results['p']< 0.0001:
            collected_results['pval_generic'] = float('{:.1e}'.format(uml_results['p']))
        else:
            collected_results['pval_generic'] = round(uml_results['p'],4)

        collected_results['pval_bayesian'] = None

        if collected_results['pval_generic']<0.01:
            send_notif=True

        if collected_results['pval_generic'] <0.1:
            uml_ontime = format_ontime_events_uml(uml_results['coincident_events'], event_mjd)
            coinc_events=[]
            for eventid in uml_ontime.keys():
                if (uml_ontime[eventid]['event_pval_generic'] < 0.1):
                    uml_ontime[eventid]['event_pval_bayesian'] = None
                    coinc_events.append(uml_ontime[eventid])
            collected_results['n_events_coincident'] = len(coinc_events)
            collected_results['coincident_events'] = coinc_events

            collected_results['most_probable_direction'] = {
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
        collected_results['pval_generic'] = None
        
        if llama_results['p_value'] < 0.0001:
            collected_results['pval_bayesian'] = float('{:.1e}'.format(llama_results['p_value']))
        else:
            collected_results['pval_bayesian'] = round(llama_results['p_value'],4)

        if collected_results['pval_bayesian']<0.01:
            send_notif=True

        if collected_results['pval_bayesian']<0.1:
            llama_ontime = format_ontime_events_llama(llama_results['single_neutrino'])
            coinc_events=[]
            for eventid in llama_ontime.keys():
                if (llama_ontime[eventid]['event_pval_bayesian'] < 0.1):
                    llama_ontime[eventid]['event_pval_generic'] = None
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
    if ('most_probable_direction' in collected_results.keys()):
        if (collected_results['n_events_coincident'] == 0):
            c = collected_results.pop('most_probable_direction')
    if ('most_probable_direction' in collected_results.keys()):
        try: 
            if (collected_results['pval_generic']>0.1):
                c = collected_results.pop('most_probable_direction')
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

        if status ==0:
            with open('/home/jthwaites/private/tokens/gw_token.txt') as f:
                my_key = f.readline()
            
            if not subthreshold:
                channels = ['#gwnu-heartbeat', '#alerts']#, '#gwnu']
            else:
                channels = ['#gwnu-heartbeat']
            for channel in channels:
                with open(os.path.join(save_location, f'{name}_collected_results.json'),'r') as fi:
                    response = requests.post('https://slack.com/api/files.upload',
                                        timeout=60,
                                        params={'token': my_key},
                                        data={'filename':'gcn.json',
                                            'title': f'GCN Notice for {name}',
                                            'channels': channel},
                                        files={'file': fi}
                                        )
                if response.ok is True:
                    logger.info("GCN posted OK to {}".format(channel))
                else:
                    logger.info("Error posting skymap to {}!".format(channel))
        
            collected_results['subthreshold'] = subthreshold
            try:
                updateGW_public(collected_results)
                logger.info('Updated webpage.')
            except:
                logger.warning('Failed to push to public webpage.')
        
        if send_notif:
            sender_script = os.path.join(os.environ.get('FAST_RESPONSE_SCRIPTS'),
                                         '../slack_posters/lvk_email_sms_notif.py')
            try: 
                subprocess.call([sender_script, '--path_to_gcn', 
                                 os.path.join(save_location, f'{name}_collected_results.json')])
                logger.info('Sent alert to ROC for p<0.01')
            except:
                logger.warning('Failed to send email/SMS notification.')

            # try:
            #     if params['Group'] == 'Burst': 
            #         merger_type = 'Burst'
            #     else:
            #         k = ['BNS','NSBH','BBH']
            #         probs = {j: float(params[j]) for j in k}
            #         merger_type = max(zip(probs.values(), probs.keys()))[1]
            # except:
            #     logger.info('Could not determine type of event')
            #     merger_type = None

            # try:
            #     subprocess.call(['/home/jthwaites/private/make_call.py', f'--type={merger_type}', '--call_anyway'])
            # except Exception as e:
            #     logger.warning('Call for p<0.01 failed.')

        else: 
            logger.info('p>0.01: no email/sms sent')
    else:
        with open(os.path.join(save_location, f'mocks/{name}_collected_results.json'),'w') as f:
            json.dump(collected_results, f, indent = 6)
        #logger.info('sending test notice')
        #status = SendTestAlert(results = collected_results)
        #logger.info('status: {}'.format(status))

        #send the notice to slack (#gw-mock-heartbeat)
        with open('/home/jthwaites/private/tokens/gw_token.txt') as f:
            my_key = f.readline()
                
        channel = '#gw-mock-heartbeat'
        with open(os.path.join(save_location, f'mocks/{name}_collected_results.json'),'r') as fi:
            response = requests.post('https://slack.com/api/files.upload',
                                    timeout=60,
                                    params={'token': my_key},
                                    data={'filename':'gcn.json',
                                          'title': f'GCN Notice for {name}',
                                          'channels': channel},
                                    files={'file': fi}
                                    )
        if response.ok is True:
            logger.info("GCN posted OK to {}".format(channel))
        else:
            logger.info("Error posting skymap to {}!".format(channel))

parser = argparse.ArgumentParser(description='Combine GW-Nu results')
parser.add_argument('--run_live', action='store_true', default=False,
                    help='Run on live GCNs')
parser.add_argument('--test_path', type=str, default=None,
                    help='path to test xml file')
parser.add_argument('--max_wait', type=float, default=35.,
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
#llama_results_location = '/home/azhang/public_html/llama/json/'
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
    
    parse_notice(record, wait_for_llama=args.wait_for_llama, heartbeat=args.heartbeat)
logger.info("done")
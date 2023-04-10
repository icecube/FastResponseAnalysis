'''Script to check for results from UML and LLAMA 
and write a gcn notice for a gw followup
written by Jessie Thwaites, March 2023
'''

import gcn, logging 
from astropy.time import Time
import time, os
import dateutil.parser
from datetime import datetime
import numpy as np
import pickle, json
from scipy.special import erfinv
import numpy as np

fra_results_location = '/home/jthwaites/FastResponse/'
llama_results_location = '/home/jthwaites/FastResponse/' #fix here, also ln 80
max_wait = 10 #minutes
save_location = '/home/jthwaites/FastResponse/' #where to save final json

@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE)

def process_gcn(payload, root):
    logger = logging.getLogger()
    logger.info("alert found, processing GCN")

    collected_results = {}
    ### GENERAL PARAMETERS ###
    #read event information
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}

    eventtime = root.find('.//ISOTime').text
    event_mjd = Time(eventtime, format='isot').mjd
    name = root.attrib['ivorn'].split('#')[1] 

    collected_results['gw_event_name'] = name.split('-')[0]
    collected_results['gw_map_num'] = int(name.split('-')[1])

    collected_results['t_merger'] = Time(eventtime, format='isot').iso
    collected_results['t_start'] = Time(event_mjd - 500./86400., format = 'mjd').iso
    collected_results['t_stop'] = Time(event_mjd + 500./86400., format = 'mjd').iso

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

        if root.attrib['role'] != 'observation':
            name=name+'_test'

        uml_results_path = os.path.join(fra_results_location, start_str + '_' + name.replace(' ', '_') \
                                   + '/' + start_str + '_' + name.replace(' ', '_')+'_results.pickle')
        uml_results_finished = os.path.exists(uml_results_path)

        #FIX HERE FOR LLAMA
        llama_results_path = os.path.join(llama_results_location,'significance_opa_lvc-i3_high.json')
        llama_results_finished = os.path.exists(llama_results_path)

        if uml_results_finished and llama_results_finished:
            results_done=True
        else:
            current_mjd = Time(datetime.utcnow(), scale='utc').mjd
            if current_mjd - start_check_mjd < max_delay:
                time.sleep(10.)
            else:
                if uml_results_finished or llama_results_finished:
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
    
    ### COLLECT RESULTS ###
    if uml_results_finished:
        with open(uml_results_path, 'rb') as f:
            results = pickle.load(f)
        
        overall_p_gen_transient = results['p']
        collected_results['sens_low'] = round(results['sens_range'][0],3)
        collected_results['sens_high'] =round(results['sens_range'][1],3)
        overall_sig_gen_transient = round(np.sqrt(2)*erfinv(1-2*overall_p_gen_transient),4)
        
        uml_coinc_nu = results['coincident_events']
        
    if llama_results_finished:
        import json
        with open(llama_results_path, 'r') as f:
            results = json.load(f)

        overall_p_bayesian = results['p_value_total']
        overall_sig_bayesian = round(np.sqrt(2)*erfinv(1-2*overall_p_bayesian),4)

        llama_coinc_nu = results['single_neutrino']

    ### COINCIDENT NEUTRINOS ###
    #if p<10% in either analysis, include coincident neutrino information
    if uml_results_finished and llama_results_finished:
        if (overall_p_gen_transient < 0.1) or (overall_p_bayesian < 0.1):
            collected_results['overall_p_bayesian'] = round(overall_p_bayesian,4)
            collected_results['overall_sig_bayesian'] = overall_sig_bayesian
            collected_results['overall_p_gen_transient'] = round(overall_p_gen_transient,4)
            collected_results['overall_sig_gen_transient'] = overall_sig_gen_transient

            ontime_events={}
            for event in uml_coinc_nu:
                ontime_events[event['event']]={
                    'dt' : round((event['time']-event_mjd)*86400., 2),
                    'ra' : round(np.rad2deg(event['ra']), 2),
                    'dec' : round(np.rad2deg(event['dec']), 2),
                    'ang_unc' : round(np.rad2deg(event['sigma']*2.145966),2),
                    'p_gen_transient': round(event['pvalue'],4)
                }
            for event in llama_coinc_nu:
                ontime_events[event['i3event']]['p_bayesian'] = round(event['p_value'],4)
            
            coinc_events=[]
            for eventid in ontime_events.keys():
                if (ontime_events[eventid]['p_gen_transient'] < 0.1) or (ontime_events[eventid]['p_bayesian'] < 0.1):
                    coinc_events.append(ontime_events[eventid])
            
            collected_results['n_events_coinc'] = len(coinc_events)
            collected_results['coinc_events'] = coinc_events

    elif uml_results_finished:
        if overall_p_gen_transient < 0.1:
            collected_results['overall_p_gen_transient'] = round(overall_p_gen_transient,4)
            collected_results['overall_sig_gen_transient'] = overall_sig_gen_transient

            coinc_events=[]
            for event in uml_coinc_nu:
                if event['pvalue'] < 0.1:
                    coinc_events.append({
                    'dt' : round((event['time']-event_mjd)*86400.,2),
                    'ra' : round(np.rad2deg(event['ra']),2),
                    'dec' : round(np.rad2deg(event['dec']),2),
                    'ang_unc' : round(np.rad2deg(event['sigma']*2.145966),2),
                    'p_gen_transient': round(event['pvalue'],4),
                    'p_bayesian': None
                    })

            collected_results['n_events_coinc'] = len(coinc_events)
            collected_results['coinc_events'] = coinc_events
    
    else:
        if overall_p_bayesian < 0.1:
            collected_results['overall_p_bayesian'] = round(overall_p_bayesian,4)
            collected_results['overall_sig_bayesian'] = overall_sig_bayesian

            coinc_events=[]
            for event in llama_coinc_nu:
                if event['p_value'] < 0.1:
                    coinc_events.append({
                    'dt' : round(event['dt'],2),
                    'ra' : round(event['ra'],2),
                    'dec' : round(event['dec'],2),
                    'ang_unc' : round(event['sigma']*2.145966,2),
                    'p_gen_transient': None,
                    'p_bayesian': round(event['p_value'],4)
                    })

            collected_results['n_events_coinc'] = len(coinc_events)         
            collected_results['coinc_events'] = coinc_events

    ### SAVE RESULTS ###
    with open(os.path.join(save_location, f'{name}_collected_results.json'),'w') as f:
        json.dump(collected_results, f, indent = 6)
    

import lxml.etree
import argparse
parser = argparse.ArgumentParser(description='Combine GW-Nu results')
parser.add_argument('--run_live', action='store_true', default=False,
                        help='Run on live GCNs')
args = parser.parse_args()

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.warning("combine_results starting, connecting to GCN")

if args.run_live:
    gcn.listen(handler=process_gcn)

else:
    logger.info('testing on mocks')
    import glob
    paths=glob.glob('/home/jthwaites/FastResponse/*/*xml')
    payload = open(paths[1], 'rb').read()
    root = lxml.etree.fromstring(payload)
    process_gcn(payload, root)

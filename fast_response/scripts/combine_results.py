'''Script to check for results from UML and LLAMA 
and write a gcn notice for a gw followup
written by Jessie Thwaites, Jan 2023
'''

import gcn
import logging
import os
from astropy.time import Time
from datetime import datetime
import time
import dateutil.parser
import numpy as np
from scipy.special import erfinv

@gcn.handlers.include_notice_types(
    #gcn.notice_types.LVC_EarlyWarning,
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE,
    gcn.notice_types.LVC_RETRACTION)

def process_gcn(payload, root):
    logger = logging.getLogger()
    logger.info("alert found, processing GCN")

    #read event information
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}

    # Read trigger time of event
    eventtime = root.find('.//ISOTime').text
    event_mjd = Time(eventtime, format='isot').mjd

    current_mjd = Time(datetime.utcnow(), scale='utc').mjd
    needed_delay = 1000./84600./2.
    current_delay = current_mjd - event_mjd

    while current_delay < needed_delay:
        logger.info("Wait {:.1f} seconds for data".format(
            (needed_delay - current_delay)*86400.)
            )
        time.sleep((needed_delay - current_delay)*86400.)
        current_mjd = Time(datetime.utcnow(), scale='utc').mjd
        current_delay = current_mjd - event_mjd

    #check for final results from uml, llama
    #uml_base_path = os.environ.get('FAST_RESPONSE_SCRIPTS') +'/'
    uml_base_path='/home/jthwaites/FastResponse/'
    if uml_base_path is None: 
        logger.info('env variable FAST_RESPONSE_SCRIPTS not set')
        logger.info('checking in /home/jthwaites/FastResponse/')
        uml_base_path='/home/jthwaites/FastResponse/'
    
    start_date = Time(dateutil.parser.parse(eventtime)).datetime
    start_str = f'{start_date.year:02d}_' \
            + f'{start_date.month:02d}_{start_date.day:02d}'
    name = root.attrib['ivorn'].split('#')[1] 

    if root.attrib['role'] != 'observation':
        name=name+'_test'

    uml_results = uml_base_path + start_str + '_' + name.replace(' ', '_') \
        + '/' + start_str + '_' + name.replace(' ', '_')+'_results.pickle'
    uml_results_finished = os.path.exists(uml_results)

    #change with correct path
    llama_base_path=uml_base_path
    llama_results = llama_base_path+'sample_llama_highsig.json'
    llama_results_finished = os.path.exists(llama_results)

    max_allowed_time = 15*60. #allowed 15 minutes to finish running
    lapsed_time = 0
    while (not (uml_results_finished and llama_results_finished) and (lapsed_time < max_allowed_time)):
        time.sleep(10)
        uml_results_finished = os.path.exists(uml_results)
        llama_results_finished = os.path.exists(llama_results)

    #if neither finished: something went wrong
    if not uml_results_finished and not llama_results_finished:
        logger.warning('error: both analyses not finished after 15 min wait')
        return
    
    #start writing gcn
    gcn_fields = [
        # name         description                                   ucd                 unit        dataType  value
        ('stream',            'Stream number',                      'meta.number',      ' ',        'float',   0),
        ('gw_event',          'GW event name from LVK',             'meta.number',      ' ',        'float',   0),
        ('gw_gcn_notice num', 'GW event notice num',                'meta.number',      ' ',        'float',   0),
        ('t_merger',          'trigger time from LVK',              'time',             'UTC',      'string',  Time(eventtime, format='isot').iso),
        ('start_time',        'start time of search (trigger-500s)','time',             'UTC',      'string',  Time(event_mjd-500./86400.,format='mjd').iso),
        ('stop_time',         'stop time of search (trigger+500s)', 'time',             'UTC',      'string',  Time(event_mjd+500./86400.,format='mjd').iso),
    ]

    #pulling saved values from both analyses
    uml_saved = {}
    if uml_results_finished:
        import pickle
        with open(uml_results, 'rb') as f:
            results = pickle.load(f)
        uml_saved['pval'] = results['p']
        uml_saved['sens_range'] = results['sens_range']

        n_coincident_events=0
        if uml_saved['pval']<0.1:
            #from FRA, calculate nsigma
            uml_saved['sigma'] = np.sqrt(2)*erfinv(1-2*results['p'])

            uml_saved['coincident_events'] =[]
            for event in results['coincident_events']:
                if event['pvalue']<=0.1:
                    ra = '{:.2f}'.format(np.rad2deg(event['ra']))
                    dec = '{:.2f}'.format(np.rad2deg(event['dec']))
                    sigma = '{:.2f}'.format(np.rad2deg(event['sigma']*2.145966))
                    dt = '{:.2f}'.format((event['time']-event_mjd)*86400.)
                    event_p = event['pvalue']
                    uml_saved['coincident_events'].append(
                        {'ra':ra,'dec':dec,'ang_unc':sigma,'dt':dt,'event_p':event_p}
                    )
                    n_coincident_events+=1
        uml_saved['n_coinc_events']=n_coincident_events

    llama_saved={}
    if llama_results_finished:
        import json
        with open(llama_results, 'r') as f:
            results = json.load(f)
        llama_saved['pval'] = results['p']
        llama_saved['sens_range'] = results['sens_range']

        n_coincident_events=0
        if llama_saved['pval']<0.1:
            #from FRA, calculate nsigma
            llama_saved['sigma'] = np.sqrt(2)*erfinv(1-2*results['p'])

            llama_saved['coincident_events'] =[]
            for event in results['coincident_events']:
                if event['pvalue']<=0.1:
                    ra = '{:.2f}'.format(np.rad2deg(event['ra']))
                    dec = '{:.2f}'.format(np.rad2deg(event['dec']))
                    sigma = '{:.2f}'.format(np.rad2deg(event['sigma']*2.145966))
                    dt = '{:.2f}'.format((event['time']-event_mjd)*86400.)
                    event_p = event['pvalue']
                    llama_saved['coincident_events'].append(
                        {'ra':ra,'dec':dec,'ang_unc':sigma,'dt':dt,'event_p':event_p}
                    )
                    n_coincident_events+=1
        llama_saved['n_coinc_events']=n_coincident_events
    
    #appending gcn fields
    if llama_results_finished and uml_results_finished:
        highest_p = uml_saved['pval'] if uml_saved['pval'] > llama_saved['pval'] else llama_saved['pval']
        highest_p_ana =uml_saved if uml_saved['pval'] > llama_saved['pval'] else uml_saved
    elif llama_results_finished:
        highest_p = llama_saved['pval']
        highest_p_ana = llama_saved
    elif uml_results_finished:
        highest_p = uml_saved['pval']
        highest_p_ana = uml_saved
    
    if highest_p > 0.1:
        logger.info('p>0.1, no significant p found')
        addl_fields = [
            # name             description                                      ucd              unit        dataType  value
            ('n_events_coinc', 'number of coincdent events',                    'meta.number',   ' ',        'float',   0),
            ('high_sens',      'highest sensitivity (E2dN/dE) at map locations','meta.number',   'GeV cm^-2','float', highest_p_ana['sens_range'][1]),
            ('low_sens',       'lowest sensitivity (E2dN/dE) at map locations', 'meta.number',   'GeV cm^-2','float', highest_p_ana['sens_range'][0]),
        ]
    else:
        logger.info('p<0.1, significant followup found')
        if llama_results_finished and uml_results_finished:
            addl_fields = [
                 # name                       description                                           ucd              unit  dataType  value
                ('n_events_coinc',            'number of coincdent events',                         'meta.number',   ' ',  'float', highest_p_ana['n_coinc_events']),
                ('overall_p_gen_transient',   'overall p for generic transient followup',           'meta.number',   ' ',  'float', uml_saved['pval']),
                ('overall_p_bayesian',        'overall p for Bayesian followup',                    'meta.number',   ' ',  'float', llama_saved['pval']),
                ('overall_sig_gen_transient', 'overall significance for generic transient followup','meta.number',   ' ',  'float', uml_saved['sigma']),
                ('overall_sig_bayesian',      'overall significance for Bayesian followup',         'meta.number',   ' ',  'float', llama_saved['sigma']),
            ]
        elif llama_results_finished:
            addl_fields = [
                 # name                       description                                           ucd              unit  dataType  value
                ('n_events_coinc',            'number of coincdent events',                         'meta.number',   ' ',  'float', highest_p_ana['n_coinc_events']),
                ('overall_p_bayesian',        'overall p for Bayesian followup',                    'meta.number',   ' ',  'float', llama_saved['pval']),
                ('overall_sig_bayesian',      'overall significance for Bayesian followup',         'meta.number',   ' ',  'float', llama_saved['sigma']),
            ]
        elif uml_results_finished:
            addl_fields = [
                 # name                       description                                           ucd              unit  dataType  value
                ('n_events_coinc',            'number of coincdent events',                         'meta.number',   ' ',  'float', highest_p_ana['n_coinc_events']),
                ('overall_p_gen_transient',   'overall p for generic transient followup',           'meta.number',   ' ',  'float', uml_saved['pval']),
                ('overall_sig_gen_transient', 'overall significance for generic transient followup','meta.number',   ' ',  'float', uml_saved['sigma']),
            ]
        
        #event_num=1
        #for event in highest_p_ana['coincident_events']:
            ##TODO: table of coinc neutrinos

        addl_fields.append(
                ('high_sens',      'highest sensitivity (E2dN/dE) at map locations','meta.number',   'GeV cm^-2','float', highest_p_ana['sens_range'][1]),
        )
        addl_fields.append(
                ('low_sens',       'lowest sensitivity (E2dN/dE) at map locations', 'meta.number',   'GeV cm^-2','float', highest_p_ana['sens_range'][0]),
        )
    
    #add together all fields to form final notice  
    what = gcn_fields+addl_fields
    
    print(what)
    #then....? TODO: next

import lxml.etree
import argparse
parser = argparse.ArgumentParser(description='Fast Response Analysis')
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
    payload = open(paths[0], 'rb').read()
    root = lxml.etree.fromstring(payload)
    process_gcn(payload, root)

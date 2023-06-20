#!/usr/bin/env python

''' Script to automatically receive GCN alerts and get LIGO skymaps 
    to run realtime neutrino follow-up

    Author: Raamis Hussain, updated by Jessie Thwaites, MJ Romfoe
    Date:   March 2023
'''

import gcn
import sys
import pickle
from dateutil.parser import parse
from dateutil.relativedelta import relativedelta

@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE)

def process_gcn(payload, root):

    AlertTime=datetime.utcnow().isoformat()
    log_file.flush()
    analysis_path = os.environ.get('FAST_RESPONSE_SCRIPTS')
    if analysis_path is None:
        try:
            import fast_response
            analysis_path = os.path.join(os.path.dirname(fast_response.__file__),'scripts/')
        except Exception as e:
            print(e)
            print('###########################################################################')
            print('CANNOT FIND ENVIRONMENT VARIABLE POINTING TO REALTIME FAST RESPONSE PACKAGE\n')
            print('You can either (1) install fast_response via pip or ')
            print('(2) put \'export FAST_RESPONSE_SCRIPTS=/path/to/fra/scripts\' in your bashrc')
            print('###########################################################################')
            log_file.flush()
            exit()

    # Read all of the VOEvent parameters from the "What" section.
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}
    name = root.attrib['ivorn'].split('#')[1]
    
    # if this is the listener for real events and it gets a mock, skip running on it
    if not mock and root.attrib['role']!='observation':
        return
    # only run on significant events
    if 'Significant' in params.keys():
        if int(params['Significant'])==0: 
            #not significant, do not run
            print(f'Found a subthreshold event {name}')
            log_file.flush()
            return
    else:
        # O3 does not have this parameter, this should only happen for testing
        print('No significance parameter found in LVK GCN.')
        log_file.flush()
    
    print('\n' +'INCOMING ALERT FOUND: ',datetime.utcnow())

    if root.attrib['role']=='observation' and not mock:
        ## Call everyone because it's a real event!
        username = pwd.getpwuid(os.getuid())[0]
        #if username == 'realtime':
        #    call_command = [os.path.join(analysis_path, 'make_call.py')]
        #else:
        #    call_command=['/cvmfs/icecube.opensciencegrid.org/users/jthwaites/make_call.py']
        call_command=['/cvmfs/icecube.opensciencegrid.org/users/jthwaites/make_call.py']
    
        call_args = ['--justin']
        for arg in call_args:
            call_command.append(arg+'=True')
        try:
            subprocess.call(call_command)
            #print('Call here.')
        except Exception as e:
            print('Call failed.')
            print(e)
            log_file.flush()
    
    # Read trigger time of event
    eventtime = root.find('.//ISOTime').text
    event_mjd = Time(eventtime, format='isot').mjd
    print(f'Alert MJD: {event_mjd}')
    print('GW merger time: %s \n' % Time(eventtime, format='isot').iso)
    log_file.flush()

    current_mjd = Time(datetime.utcnow(), scale='utc').mjd
    needed_delay = 1000./84600./2.
    current_delay = current_mjd - event_mjd

    # We need to make sure all the data has been collected before we can run.
    # Check to see if we need to wait for the +500 sec of data to arrive
    FiveHundred_delay = (needed_delay - current_delay)*86400.

    while current_delay < needed_delay:
        print("Need to wait another {:.1f} seconds before running".format(
            (needed_delay - current_delay)*86400.)
            )
        log_file.flush()
        time.sleep((needed_delay - current_delay)*86400.)
        current_mjd = Time(datetime.utcnow(), scale='utc').mjd
        current_delay = current_mjd - event_mjd

    skymap = params['skymap_fits']

    # Skymap distributed is a different format than previous, but previous is available.
    # Download correct format from GraceDB
    if 'multiorder' in skymap:
        time.sleep(6.) #if we don't wait, the old format isn't uploaded
        try:
            #HERE: added .fits but not impl yet. need to keep baystar/bilby naming?
            bayestar_map=skymap.replace('multiorder.','').split(',')
            if len(bayestar_map)==1:
                suffix='.gz'
            else: 
                suffix = '.gz,'+ bayestar_map[1]
            new_map = bayestar_map[0]+suffix
            map_type= bayestar_map[0].split('/')[-1]
            
            wget.download(new_map, out=os.environ.get('FAST_RESPONSE_OUTPUT')+f'skymaps/{name}_{map_type}{suffix}')
            skymap=os.environ.get('FAST_RESPONSE_OUTPUT')+f'skymaps/{name}_{map_type}{suffix}'
        except:
            print('Failed to download skymap in correct format! \nDownload skymap and then re-run script with')
            print(f'args:  --time {event_mjd} --name {name} --skymap PATH_TO_SKYMAP')
            return
            #TODO: add make_call to me here if this fails

    if root.attrib['role'] != 'observation':
        name=name+'_test'
        print('Running on scrambled data')
        log_file.flush()
    command = os.path.join(analysis_path, 'run_gw_followup.py')

    print('Running {}'.format(command))
    log_file.flush()

    subprocess.call([command, '--skymap={}'.format(skymap), 
        '--time={}'.format(str(event_mjd)), 
        '--name={}'.format(name)]
        #'--allow_neg_ts=True']
        )
    
    output = os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'),eventtime[0:10].replace('-','_')+'_'+name)
    #update webpages
    webpage_update = os.path.join(analysis_path,'document.py')
    if not mock and root.attrib['role'] == 'observation':
        try:
            subprocess.call([webpage_update,  '--gw', f'--path={output}'])

            wp_link = 'https://user-web.icecube.wisc.edu/~jthwaites/FastResponse/gw-webpage/output/{}.html'.format(eventtime[0:10].replace('-','_')+'_'+name)
            slack_message = {
                "icon_emoji": ":gw:",
                "text" : "UML GW analysis finished running for event {}: <{}|link>.".format(name, wp_link),
            } 
            with open('../slack_posters/internal_alert_slackbot.txt') as f:
                channel = f.readline().rstrip('\n')
                webhook = f.readline().rstrip('\n')
                bot_name = f.readline().rstrip('\n')

            bot = slackbot(channel, bot_name, webhook)
            bot.send_message(slack_message['text'],slack_message['icon_emoji'])
            
        except Exception as e:
            print('Failed to push to (private) webpage.')
            print(e)
            log_file.flush()

    endtime=datetime.utcnow().isoformat()
    alert_mjd = Time(AlertTime, format='isot').mjd
    end_mjd = Time(endtime, format='isot').mjd

    # Calculate latency benchmarks and save in a pickle file
    Ligo_late_sec = (alert_mjd - event_mjd)*86400
    Ice_late_sec = (end_mjd - alert_mjd)*86400
    Total_late_sec = Ligo_late_sec + Ice_late_sec

    gw_latency = {'Trigger_Time': event_mjd, 'GCN_Alert': alert_mjd, 'End_Time': end_mjd,
                    'Ligo_Latency': Ligo_late_sec, 'IceCube_Latency': Ice_late_sec, 'Total_Latency': Total_late_sec,
                    'We_had_to_wait:': FiveHundred_delay}
    
    save_dir = 'latency_o4' if root.attrib['role']=='observation' else 'PickledMocks'

    #check for directory to save pickle files and create if needed 
    if not os.path.exists(os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'),save_dir)):
        os.mkdir(os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'),save_dir))

    with open(os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'), f'{save_dir}/gw_latency_dict_{name}.pickle'), 'wb') as file:
        pickle.dump(gw_latency, file, protocol=pickle.HIGHEST_PROTOCOL)
    
    #save xml and skymap, for later
    et = lxml.etree.ElementTree(root)
    et.write(os.path.join(output, '{}-{}-{}.xml'.format(params['GraceID'], 
                        params['Pkt_Ser_Num'], params['AlertType'])), pretty_print=True)
    
    #skymap_filename=skymap.split('/')[-1]
    #subprocess.call(['wget', skymap])
    #subprocess.call(['mv', skymap_filename, output])

    if root.attrib['role'] != 'observation':
        # Move mocks to a seperate folder to avoid swamping FRA output folder
        subprocess.call(['mv',output, '/data/user/jthwaites/o4-mocks/'])
        output = '/data/user/jthwaites/o4-mocks/' + eventtime[0:10].replace('-','_')+'_'+name
    
    print('Output directory: ',output)
    log_file.flush()

if __name__ == '__main__':
    import os, subprocess, pwd
    import healpy as hp
    import numpy as np
    import lxml.etree
    import argparse
    import time
    from astropy.time import Time
    from datetime import datetime
    from fast_response.slack_posters.slack import slackbot
    import wget

    output_path = '/home/jthwaites/public_html/FastResponse/'
    #output_path=os.environ.get('FAST_RESPONSE_OUTPUT')
    #if output_path==None:
    #    output_path=os.getcwd()

    parser = argparse.ArgumentParser(description='FRA GW followup')
    parser.add_argument('--run_live', action='store_true', default=False,
                        help='Run on live GCNs')
    parser.add_argument('--heartbeat', action = 'store_true', default=False,
                        help='Run the listener as a heartbeat, running on mock LVK events only (default=False)')
    parser.add_argument('--log_path', default=output_path, type=str,
                        help='Redirect output to a log file with this path')
    parser.add_argument('--test_path', default='S191216ap_update.xml', type=str,
                        help='Skymap for use in testing listener')
    parser.add_argument('--test_o3', default=False, action='store_true',
                        help='bool to decide if we should run an already unblinded skymap with unblinded data')
    args = parser.parse_args()

    if args.heartbeat:
        logfile=os.path.join(args.log_path,'mock_log.log')
    else:
        logfile=os.path.join(args.log_path,'log.log')

    print(f'Logging to file: {logfile}')
    original_stdout=sys.stdout
    log_file = open(logfile, "a+")
    sys.stdout=log_file
    sys.stderr=log_file

    if args.run_live:
        print("Listening for GCNs . . . ")
        log_file.flush()

        mock=args.heartbeat
        print('Starting heartbeat listener') if mock else print('Running on REAL events only')
        log_file.flush()
        
        gcn.listen(handler=process_gcn)

    else: 
        print("Offline testing . . . ")
        log_file.flush()
        
        ### FOR OFFLINE TESTING
        try:
            import fast_response
            #sample_skymap_path='/data/user/jthwaites/o3-gw-skymaps/'
            sample_skymap_path=os.path.join(os.path.dirname(fast_response.__file__),'sample_skymaps/')
        except Exception as e:
            print(e)
            sample_skymap_path='/data/user/jthwaites/o3-gw-skymaps/'
        
        payload = open(os.path.join(sample_skymap_path,args.test_path), 'rb').read()
        root = lxml.etree.fromstring(payload) 

        mock=args.heartbeat
        #test runs on scrambles, observation runs on unblinded data
        if not args.test_o3:
            root.attrib['role']='test'
            mock=True

        process_gcn(payload, root)

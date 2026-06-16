#!/usr/bin/env python

''' Script to automatically receive GCN alerts and get LIGO skymaps 
    to run realtime neutrino follow-up

    Author: Raamis Hussain, Jessie Thwaites, MJ Romfoe
    Updated Date: June 2026
'''

#import gcn
import logging
from gcn_kafka import Consumer
import sys, pickle, os, subprocess, pwd
from dateutil.parser import parse
from dateutil.relativedelta import relativedelta
import healpy as hp
import numpy as np
import lxml.etree
import argparse, time, wget
from astropy.time import Time
from datetime import datetime
from fast_response.slack_posters.slack import slackbot

print("Connecting to GCN as Consumer")

with open('/home/jthwaites/private/tokens/kafka_token.txt') as f:
    client_id = f.readline().rstrip('\n')
    client_secret = f.readline().rstrip('\n')

consumer = Consumer(client_id=client_id,
                    client_secret=client_secret,
                    domain='gcn.nasa.gov',
                    #config={'max.poll.interval.ms':1800000},
                   )

consumer.subscribe(['gcn.classic.voevent.LVC_EARLY_WARNING',
                    'gcn.classic.voevent.LVC_INITIAL',
                    'gcn.classic.voevent.LVC_PRELIMINARY',
                    'gcn.classic.voevent.LVC_RETRACTION',
                    #'gcn.classic.voevent.LVC_TEST',
                    'gcn.classic.voevent.LVC_UPDATE'])

# make a new logging.FileHandler that can flush as we go
class LogFileWriter(logging.FileHandler):
    '''Make a new logging handler that uses a file
    and flushes all messges to the file as they are written
    '''
    def emit(self, record):
        super().emit(record)
        self.flush()

def process_gcn(record): #payload, root):
    AlertTime=datetime.utcnow().isoformat()
    analysis_path = os.environ.get('FAST_RESPONSE_SCRIPTS')

    if analysis_path is None:
        try:
            import fast_response
            analysis_path = os.path.join(os.path.dirname(fast_response.__file__),'scripts/')
        except Exception as e:
            logger.error('Error finding FRA package!!')
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
    name = record.attrib['ivorn'].split('#')[1]
    
    # only run on significant events
    if 'Significant' in params.keys():
        if int(params['Significant'])==0: 
            #not significant, do not run
            logger.warning(f'Found a subthreshold event {name}')
            record.attrib['role']='test'
    else:
        # O3 does not have this parameter, this should only happen for testing
        logger.warning('No significance parameter found in LVK GCN.')
    # if this is the listener for real events and it gets a mock (or low signficance), skip it
    if not mock and record.attrib['role']!='observation':
        return
    
    logger.warning('\n' +'INCOMING ALERT FOUND: ',datetime.utcnow())

    #get type of event (burst, bbh, nsbh, bns)
    try:
        if params['Group'] == 'Burst': 
            merger_type = 'Burst'
        elif params['Search'] == 'SSM':
            merger_type='SSM'
        else:
            k = ['BNS','NSBH','BBH']
            probs = {j: float(params[j]) for j in k}
            merger_type = max(zip(probs.values(), probs.keys()))[1]
    except:
        logger.warning('Could not determine type of event')
        merger_type = None
    
    if record.attrib['role']=='observation' and not mock:
        ## Call everyone because it's a real event!
        call_command=['/home/jthwaites/private/make_call.py', f'--name={name}']
    
        call_args = ['--justin']
        for arg in call_args:
            call_command.append(arg+'=True')
        if merger_type is not None:
            call_command.append(f'--type={merger_type}')
            
        try:
            subprocess.call(call_command)
        except Exception as e:
            logger.error('Call failed!')
            logger.error(e)
            
    # want heartbeat listener not to run on real events, otherwise it overwrites the main listener output
    if mock and record.attrib['role']=='observation':
        logger.info('Listener in heartbeat mode found real event. Skipping...')
        return
    
    # Read trigger time of event
    eventtime = record.find('.//ISOTime').text
    event_mjd = Time(eventtime, format='isot').mjd
    logger.info(f'Alert MJD: {event_mjd}')
    logger.info('GW merger time: %s \n' % Time(eventtime, format='isot').iso)

    current_mjd = Time(datetime.utcnow(), scale='utc').mjd
    needed_delay = 1000./84600./2.
    current_delay = current_mjd - event_mjd

    # We need to make sure all the data has been collected before we can run.
    # Check to see if we need to wait for the +500 sec of data to arrive
    FiveHundred_delay = (needed_delay - current_delay)*86400.

    while current_delay < needed_delay:
        logger.info("Need to wait another {:.1f} seconds before running".format(
            (needed_delay - current_delay)*86400.)
            )
        time.sleep((needed_delay - current_delay)*86400.)
        current_mjd = Time(datetime.utcnow(), scale='utc').mjd
        current_delay = current_mjd - event_mjd

    skymap = params['skymap_fits']

    # Multiorder Coverage (MOC) map links are distributed over the GCNs. 
    # Download flattened (normal healpy) map from GraceDB
    if 'multiorder' in skymap:
        time.sleep(6.) #if we don't wait, the old format isn't uploaded
        try: 
            #Try to get flat-res healpy map
            flat_map=skymap.replace('multiorder.','').split(',')
            if len(flat_map)==1:
                suffix='.gz'
            else: 
                suffix = '.gz,'+ flat_map[1]
            new_map = flat_map[0]+suffix
            map_type= flat_map[0].split('/')[-1]

            wget.download(new_map, out=os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'),f'skymaps/{name}_{map_type}{suffix}'))
            skymap=os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'),f'skymaps/{name}_{map_type}{suffix}')
        except:
            logger.warning('Failed to download flat-resolution skymap. Trying to convert MOC map')

            try:
                filename=skymap.split('/')[-1]
                new_output = os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'),f'skymaps/{name}_{filename}')
                wget.download(skymap, out=new_output)
                subprocess.call([os.path.join(analysis_path, 'convert_moc_to_healpix.py'),
                                '--skymap', new_output])
                if os.path.exists(new_output.replace('multiorder','converted')):
                    skymap = new_output.replace('multiorder','converted')
                    logger.info('Successfully converted map: {}'.format(skymap))
                else:
                    raise Exception('Failed to convert map.')
            except:
                logger.error('Failed to get skymap in correct format! \nDownload skymap and then re-run script with' +\
                            f'args:  --time {event_mjd} --name {name} --skymap PATH_TO_SKYMAP')
                return

    if record.attrib['role'] != 'observation':
        name=name+'_test'
        logger.info('Running on scrambled data')
    command = os.path.join(analysis_path, 'run_gw_followup.py')

    logger.info('Running {}'.format(command))
    #### FOR NOW: testing
    return

    subprocess.call([command, '--skymap={}'.format(skymap), 
        '--time={}'.format(str(event_mjd)), 
        '--name={}'.format(name)]
        #'--allow_neg_ts=True']
        )
    
    analysis_start = Time(event_mjd - 500./86400., format='mjd').iso
    output = os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'),
                          analysis_start[0:10].replace('-','_')+'_'+name)
    #update webpages
    webpage_update = os.path.join(analysis_path,'document.py')
    if not mock and record.attrib['role'] == 'observation':
        try:
            subprocess.call([webpage_update,  '--gw', f'--path={output}'])

            wp_link = 'https://user-web.icecube.wisc.edu/~jthwaites/FastResponse/gw-webpage/output/'+\
                      '{}.html'.format(analysis_start[0:10].replace('-','_')+'_'+name)
            slack_message = "UML GW analysis finished running for event {}: <{}|link>.".format(name, wp_link) 

            for channel in ['#fra-shifting','#gwnu-heartbeat']:
                bot = slackbot(channel)
                bot.post_short_msg(slack_message)
            
        except Exception as e:
            logger.error('Failed to push to (private) webpage.')
            logger.error(e)

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
    
    save_dir = 'latency_o4' if record.attrib['role']=='observation' else 'PickledMocks'

    #check for directory to save pickle files and create if needed 
    if not os.path.exists(os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'),save_dir)):
        os.mkdir(os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'),save_dir))

    #check to see if results file was created, before saving pickle file
    #essentially, check if the analysis finished sucessfully
    if os.path.exists(os.path.join(output, '{}_{}_results.pickle'.format(eventtime[0:10].replace('-','_'), name))):
        with open(os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'), f'{save_dir}/gw_latency_dict_{name}.pickle'), 'wb') as file:
            pickle.dump(gw_latency, file, protocol=pickle.HIGHEST_PROTOCOL)
    
    #save xml and skymap, for later
    et = lxml.etree.ElementTree(record)
    et.write(os.path.join(output, '{}-{}-{}.xml'.format(params['GraceID'], 
                        params['Pkt_Ser_Num'], params['AlertType'])), pretty_print=True)

    if record.attrib['role'] != 'observation':
        # Move mocks to a seperate folder to avoid swamping FRA output folder
        subprocess.call(['mv',output, '/data/user/jthwaites/o4-mocks/'])
        output = '/data/user/jthwaites/o4-mocks/' + eventtime[0:10].replace('-','_')+'_'+name
    
    logger.info('Output directory: ',output)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='FRA GW followup')
    parser.add_argument('--run_live', action='store_true', default=False,
                        help='Run on live GCNs')
    parser.add_argument('--heartbeat', action = 'store_true', default=False,
                        help='Run the listener as a heartbeat, running on mock LVK events only (default=False)')
    parser.add_argument('--log_path', default='/home/jthwaites/public_html/FastResponse/', type=str,
                        help='Redirect output to a log file with this path. Note: this is only used when running live')
    parser.add_argument('--test_path', default='S191216ap_update.xml', type=str,
                        help='Skymap for use in testing listener')
    parser.add_argument('--test_o3', default=False, action='store_true',
                        help='bool to decide if we should run an already unblinded skymap with unblinded data')
    args = parser.parse_args()

    if args.heartbeat:
        logfile=os.path.join(args.log_path,'mock_log.log')
    else:
        logfile=os.path.join(args.log_path,'log.log')

    if args.run_live:
        print(f'Logging to file: {logfile}')
         
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)
        # adds the logfile as an additional logger. this will also log to stout
        logger.addHandler(LogFileWriter(logfile, mode='a+'))
        logger.warning("Listening for GCNs . . . ")

        mock=args.heartbeat
        logger.info('Starting heartbeat listener') if mock else logger.info('Running on REAL events only')
        
        try:
            while True:
                for message in consumer.consume(timeout=1):
                    if message.error():
                        logger.warning(message.error())
                        continue
                    value = message.value()
                    logger.warning('Found GCN on topic {}'.format(message.topic()))
                    notice = lxml.etree.fromstring(value.decode('utf-8').encode('ascii'))
                    parse_notice(notice)
        except KeyboardInterrupt:
            # make sure the logfile gets shutdown correctly and file closed
            logger.shutdown()

    else:
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)
        logger.warning("Offline testing . . . ")
       
        ### FOR OFFLINE TESTING
        try:
            import fast_response
            #sample_skymap_path='/data/user/jthwaites/o3-gw-skymaps/'
            sample_skymap_path=os.path.join(os.path.dirname(fast_response.__file__),'sample_skymaps/')
        except Exception as e:
            logger.error('Failed to find sample skymap paths!')
            logger.error(e)
            sample_skymap_path='/data/user/jthwaites/o3-gw-skymaps/'
        
        #payload = open(os.path.join(sample_skymap_path,args.test_path), 'rb').read()
        payload = open(args.test_path,'rb').read()
        record = lxml.etree.fromstring(payload) 

        mock=args.heartbeat
        #test runs on scrambles, observation runs on unblinded data
        if not args.test_o3:
            record.attrib['role']='test'
            mock=True

        process_gcn(record)

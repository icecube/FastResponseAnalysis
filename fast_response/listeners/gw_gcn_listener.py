''' Script to automatically receive GCN alerts and get LIGO skymaps 
    to run realtime neutrino follow-up

    Author: Raamis Hussain, updated by Jessie Thwaites, MJ Romfoe
    Date:   October 2022
'''

import gcn
import sys

@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE)

def process_gcn(payload, root):

    print('INCOMING ALERT FOUND: ',datetime.utcnow())
    AlertTime=datetime.utcnow().isoformat()
    log_file.flush()
    analysis_path = os.environ.get('FAST_RESPONSE_SCRIPTS')
    if analysis_path is None:
        try:
            import fast_response
            analysis_path = os.path.dirname(fast_response.__file__) + '/scripts/'
        except Exception as e:
            print(e)
            print('###########################################################################')
            print('CANNOT FIND ENVIRONMENT VARIABLE POINTING TO REALTIME FAST RESPONSE PACKAGE\n')
            print('You can either (1) install fast_response via pip or ')
            print('(2) put \'export FAST_RESPONSE_SCRIPTS=/path/to/fra/scripts\' in your bashrc')
            print('###########################################################################')
            log_file.flush()
            exit()

    # Respond only to 'test' events.
    # VERY IMPORTANT! Replace 'test' with
    # 'observation' to repond to real events.
    #if root.attrib['role']=='observation':
        ## Call everyone because it's a real event
    #    call_command = analysis_path + '/make_call.py'
    #    call_args = ['--raamis','--justin','--alex']
    #    subprocess.call([call_command,call_args[0],'True',call_args[1],'True',call_args[2],'True'])

    # Read all of the VOEvent parameters from the "What" section.
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}

    # Read trigger time of event
    eventtime = root.find('.//ISOTime').text
    event_mjd = Time(eventtime, format='isot').mjd

    print('GW trigger time: %s \n' % Time(eventtime, format='isot').iso)
    log_file.flush()

    current_mjd = Time(datetime.utcnow(), scale='utc').mjd
    needed_delay = 1000./84600./2.
    current_delay = current_mjd - event_mjd
    while current_delay < needed_delay:
        print("Need to wait another {:.1f} seconds before running".format(
            (needed_delay - current_delay)*86400.)
            )
        log_file.flush()
        time.sleep((needed_delay - current_delay)*86400.)
        current_mjd = Time(datetime.utcnow(), scale='utc').mjd
        current_delay = current_mjd - event_mjd

    skymap = params['skymap_fits']
    name = root.attrib['ivorn'].split('#')[1]

    if root.attrib['role'] != 'observation':
        name=name+'_test'
        print('Running on scrambled data')
        log_file.flush()
    command = analysis_path + 'run_gw_followup.py'

    print('Running {}'.format(command))
    log_file.flush()

    subprocess.call([command, '--skymap={}'.format(skymap), 
        '--time={}'.format(str(event_mjd)), 
        '--name={}'.format(name)]
        #'--allow_neg_ts=True']
        )
    endtime=datetime.utcnow().isoformat()

    from dateutil.parser import parse
    from dateutil.relativedelta import relativedelta

#Creates txt file for latency evaluation
    file_object = open('MilestoneTimes.txt', "a+")
    file_object.write("GCN Alert" +repr(AlertTime) +'\n' +"Trigger Time" +repr(eventtime) +'\n' +"End Time" +repr(endtime))
    file_object.close()

    file_object = open('GWLatency.txt', "a+")
    time_1 = parse(eventtime)
    time_2 = parse(AlertTime)
    time_3 = parse(endtime)

    delta_Ligo = relativedelta(time_2, time_1)
    delta_Ice = relativedelta(time_3, time_2)
    delta_total = relativedelta(time_3, time_1)

    file_object.write("Ligo Latency:" +repr(delta_Ligo) +'\n' +"IceCube Latency:" +repr(delta_Ice) +'\n' + "Total Latency:" +repr(delta_total))
    file_object.close()


    for directory in os.listdir(analysis_path+'../../output'):
        if name in directory: 
            skymap_filename=skymap.split('/')[-1]
            if 'MS22' in name:
                #import glob
                et = lxml.etree.ElementTree(root)
                et.write(analysis_path+'../../output/'+directory+'/{}-{}-{}.xml'.format(params['GraceID'], 
                         params['Pkt_Ser_Num'], params['AlertType']), pretty_print=True)
                subprocess.call(['wget', skymap])
                subprocess.call(['mv', skymap_filename, analysis_path+'../../output/'+directory])
                subprocess.call(['mv',analysis_path+'../../output/'+directory, '/data/user/jthwaites/o4-mocks/'])
                print('Output directory: ','/data/user/jthwaites/o4-mocks/'+directory)
                log_file.flush()
            else:
                print('Output directory: ',analysis_path+'../../output/'+directory)
                log_file.flush()
            break

if __name__ == '__main__':
    import os, subprocess
    import healpy as hp
    import numpy as np
    import lxml.etree
    import argparse
    import time
    from astropy.time import Time
    from datetime import datetime

    output_path=os.environ.get('FAST_RESPONSE_OUTPUT')
    if output_path==None:
        output_path=os.getcwd()

    parser = argparse.ArgumentParser(description='FRA GW followup')
    parser.add_argument('--run_live', action='store_true', default=False,
                        help='Run on live GCNs')
    parser.add_argument('--log_path', default=output_path, type=str,
                        help='Redirect output to a log file with this path')
    args = parser.parse_args()
    logfile=args.log_path
    original_stdout=sys.stdout
    log_file = open(logfile, "a+")
    sys.stdout=log_file
    sys.stderr=log_file

    print("Wowza!!!")
    log_file.flush()
          

    if args.run_live:
        print("Listening for GCNs . . . ")
        log_file.flush()
        gcn.listen(handler=process_gcn)
    else: 
        print("Listening for GCNs . . . ")
        log_file.flush()
        ### FOR OFFLINE TESTING
        try:
            import fast_response
            sample_skymap_path='/data/user/jthwaites/o3-gw-skymaps/'
            #sample_skymap_path=os.path.dirname(fast_response.__file__) +'/sample_skymaps/'
        except Exception as e:
            sample_skymap_path='/data/user/jthwaites/o3-gw-skymaps/'
        
        payload = open(sample_skymap_path + 'S190728q-5-Update.xml', 'rb').read()
        root = lxml.etree.fromstring(payload) 

        #test runs on scrambles, observation runs on unblinded data
        #root.attrib['role']='test'
        process_gcn(payload, root)

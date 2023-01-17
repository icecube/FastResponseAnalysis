#!/usr/bin/env python

''' Script to run a 2 week follow-up for a GW event containing a NS
    adapted from gw_gcn_listener.py

    Author: Jessie Thwaites
    Date:   January 2023
'''
#import gcn

def process_gcn(payload, root):
    print('Running 2 week followup for GW with NS')
    print('Starting: ',datetime.utcnow())
    
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
            exit()

    # Read all of the VOEvent parameters from the "What" section.
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}

    # Read trigger time of event
    if root.attrib['role'] == 'observation':
        eventtime = root.find('.//ISOTime').text
        event_mjd = Time(eventtime, format='isot').mjd
    else:
        #if testing, want to query livestream rather than load archival, so use recent time
        eventtime = '2022-12-15T21:51:25.506'
        event_mjd = Time(eventtime, format='isot').mjd

    print('GW merger time: %s \n' % Time(eventtime, format='isot').iso)

    skymap = params['skymap_fits']
    name = root.attrib['ivorn'].split('#')[1]

    # Skymap distributed is a different format than previous, but previous is available.
    # Download correct format from GraceDB
    if 'multiorder' in skymap:
        try:
            bayestar_map=skymap.replace('multiorder.','').split(',')
            if len(bayestar_map)==1:
                new_map=bayestar_map[0]+'.gz'
            else: 
                new_map=bayestar_map[0]+'.gz,'+ bayestar_map[1]
            import requests
            ret = requests.head(new_map)
            assert ret.status_code == 200
            skymap=new_map
        except:
            print('Failed to download skymap in correct format')
            #TODO: add make_call to me here if this fails

    name=name+'_2week'
    if root.attrib['role'] != 'observation':
        name=name+'_test'
        print('Running on scrambled data')
    command = analysis_path + 'run_gw_followup.py'

    print('Running {}'.format(command))

    subprocess.call([command, '--skymap={}'.format(skymap), 
        '--time={}'.format(str(event_mjd)), 
        '--name={}'.format(name),
        '--tw=1218240'] #[-0.1, +14]day
        #'--allow_neg_ts=True']
        )

    if os.environ.get('FAST_RESPONSE_OUTPUT') is not None:
        print('Finished \n Results saved to '+os.environ.get('FAST_RESPONSE_OUTPUT'))
    else: 
        print('Finished')

if __name__ == '__main__':
    import os, subprocess
    import numpy as np
    import lxml.etree
    import argparse
    from astropy.time import Time
    from datetime import datetime

    output_path=os.environ.get('FAST_RESPONSE_OUTPUT')
    if output_path==None:
        output_path=os.getcwd()

    parser = argparse.ArgumentParser(description='FRA GW followup')
    parser.add_argument('--path', default='../sample_skymaps/S191216ap_update.xml', type=str,
                        help='path to xml skymap (default S191216ap test)')
    #Test bool default is true for now, but will be changed to False during O4
    parser.add_argument('--test', default=True, action='store_true',
                        help='bool to run on scrambled data')
    args = parser.parse_args()

    payload = open(args.path, 'rb').read()
    root = lxml.etree.fromstring(payload) 

    if args.test:
        print("Offline testing . . . ")
        root.attrib['role']='test'
    else: 
        print('Running on unscrambled data')

    process_gcn(payload, root)
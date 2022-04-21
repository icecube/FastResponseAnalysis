''' Script to automatically receive GCN alerts and get LIGO skymaps 
    to run realtime neutrino follow-up

    Author: Raamis Hussain, updated by Jessie Thwaites
    Date:   April 2022
'''

import gcn

@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE)
def process_gcn(payload, root):

    print('INCOMING ALERT FOUND')
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

    current_mjd = Time(datetime.utcnow(), scale='utc').mjd
    needed_delay = 1000./84600./2.
    current_delay = current_mjd - event_mjd
    while current_delay < needed_delay:
        print("Need to wait another {:.1f} seconds before running".format(
            (needed_delay - current_delay)*86400.)
            )
        time.sleep((needed_delay - current_delay)*86400.)
        current_mjd = Time(datetime.utcnow(), scale='utc').mjd
        current_delay = current_mjd - event_mjd

    skymap = params['skymap_fits']
    name = root.attrib['ivorn'].split('#')[1]
    name= name.split('-')[0] #remove? this removes the map identifiers from the name

    if root.attrib['role'] != 'observation':
        name=name+'_test'
    command = analysis_path + 'run_gw_followup.py'

    print('Running {}'.format(command))
    subprocess.call([command, '--skymap={}'.format(skymap), 
        '--time={}'.format(str(event_mjd)), 
        '--name={}'.format(name)]
        )


if __name__ == '__main__':
    import os, subprocess
    import healpy as hp
    import numpy as np
    import lxml.etree
    import argparse
    import time
    from astropy.time import Time
    from datetime import datetime

    parser = argparse.ArgumentParser(description='FRA GW followup')
    parser.add_argument('--run_live', action='store_true', default=False,
                        help='Run on live GCNs')
    args = parser.parse_args()

    if args.run_live:
        print("Listening for GCNs . . . ")
        gcn.listen(handler=process_gcn)
    else: 
        ### FOR OFFLINE TESTING
        try:
            import fast_response
            sample_skymap_path=os.path.dirname(fast_response.__file__) +'/sample_skymaps/'
        except Exception as e:
            sample_skymap_path='/data/user/apizzuto/fast_response_skylab/' \
                + 'fast-response/fast_response/sample_skymaps/'
        
        payload = open(sample_skymap_path \
            + 'S191216ap_update.xml', 'rb').read()
        root = lxml.etree.fromstring(payload)

        #test runs on scrambles, observation runs on unblinded data
        root.attrib['role']='test'

        process_gcn(payload, root)
        ### For more test skymaps:
        # https://gracedb.ligo.org/api/superevents/MS190403g/files/bayestar.fits.gz
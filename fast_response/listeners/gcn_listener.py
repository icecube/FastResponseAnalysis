''' Script to automatically receive GCN notices for IceCube
    alert events and run followup accordingly

    Author: Alex Pizzuto
    Date:   July 2020
'''

import gcn
@gcn.handlers.include_notice_types(
        gcn.notice_types.ICECUBE_ASTROTRACK_GOLD,
        gcn.notice_types.ICECUBE_ASTROTRACK_BRONZE,
        gcn.notice_types.ICECUBE_CASCADE)

def process_gcn(payload, root):
    print("INCOMING ALERT")
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

    stream = params['Stream']
    if stream == '26':
        print("Detected cascade type alert, running cascade followup. . . ")
        alert_type='cascade'
        skymap = params['skymap_fits']
    else:
        print("Found track type alert, running track followup. . . ")
        alert_type='track'

    event_id = params['event_id']
    run_id = params['run_id']
    eventtime = root.find('.//ISOTime').text
    event_mjd = Time(eventtime, format='isot').mjd

    if alert_type == 'cascade':
        command = analysis_path + 'run_cascade_followup.py'
    else:
        command = analysis_path + 'run_track_followup.py'

    current_mjd = Time(datetime.utcnow(), scale='utc').mjd
    needed_delay = 1.
    current_delay = current_mjd - event_mjd
    while current_delay < needed_delay:
        print("Need to wait another {:.1f} seconds before running".format(
            (needed_delay - current_delay)*86400.)
            )
        time.sleep((needed_delay - current_delay)*86400.)
        current_mjd = Time(datetime.utcnow(), scale='utc').mjd
        current_delay = current_mjd - event_mjd

    if alert_type == 'track':
        base_skymap_path = '/home/followup/output_plots/'
        skymap_f = glob(base_skymap_path \
            + f'run{int(run_id):08d}.evt{int(event_id):012d}.*.fits.gz')
        if len(skymap_f) == 0:
            print("COULD NOT FIND THE SKYMAP FILE FOR V2 TRACK ALERT EVENT")
            return
        elif len(skymap_f) == 1:
            skymap = skymap_f[0]
        else:
            print("TOO MANY OPTIONS FOR THE SKYMAP FILE FOR V2 TRACK ALERT EVENT")
            return

    subprocess.call([command, '--skymap={}'.format(skymap), 
        '--time={}'.format(str(event_mjd)), 
        '--alert_id={}'.format(run_id+':'+event_id)]
        )

if __name__ == '__main__':
    import os, subprocess
    import healpy as hp
    import numpy as np
    import lxml.etree
    import argparse
    from astropy.time import Time
    from datetime import datetime
    import time
    from glob import glob

    parser = argparse.ArgumentParser(description='Fast Response Analysis')
    parser.add_argument('--run_live', action='store_true', default=False,
                        help='Run on live GCNs')
    parser.add_argument('--test_cascade', default=False, action='store_true',
                        help='When testing, raise to run a cascade, else track')
    args = parser.parse_args()

    if args.run_live:
        print("Listening for GCNs . . . ")
        gcn.listen(handler=process_gcn)
    elif not args.test_cascade:
        print("Running on sample track . . . ")
        payload = open('/data/user/apizzuto/fast_response_skylab/' \
            + 'fast-response/fast_response/sample_skymaps/' \
            + 'sample_astrotrack_alert_2021.xml', 'rb').read()
        root = lxml.etree.fromstring(payload)
        process_gcn(payload, root)
    else:
        print("Running on sample cascade . . . ")
        payload = open('/data/user/apizzuto/fast_response_skylab/' \
            + 'fast-response/fast_response/sample_skymaps/' \
            + 'sample_cascade.txt', 'rb').read()
        root = lxml.etree.fromstring(payload)
        process_gcn(payload, root)

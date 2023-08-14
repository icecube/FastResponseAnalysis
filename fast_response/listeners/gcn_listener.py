''' Script to automatically receive GCN notices for IceCube
    alert events and run followup accordingly

    Author: Alex Pizzuto
    Date:   July 2020
'''

from itertools import count
import gcn
@gcn.handlers.include_notice_types(
        gcn.notice_types.ICECUBE_ASTROTRACK_GOLD,
        gcn.notice_types.ICECUBE_ASTROTRACK_BRONZE,
        gcn.notice_types.ICECUBE_CASCADE)

def process_gcn(payload, root):
    print("INCOMING ALERT: ",datetime.utcnow())

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
    eventtime = root.find('.//ISOTime').text
    if stream == '26':
        print("Detected cascade type alert, running cascade followup. . . ")
        alert_type='cascade'
        event_name='IceCube-Cascade_{}{}{}'.format(eventtime[2:4],eventtime[5:7],eventtime[8:10])

        skymap = params['skymap_fits']
    else:
        print("Found track type alert, running track followup. . . ")
        alert_type='track'
        event_name='IceCube-{}{}{}'.format(eventtime[2:4],eventtime[5:7],eventtime[8:10]) 

        # IceCube sends 2: a notice and a revision, only want to run once
        if int(params['Rev']) !=0:
            return        

    event_id = params['event_id']
    run_id = params['run_id']
    event_mjd = Time(eventtime, format='isot').mjd
    try:
        bot.send_message(f'Listener found {alert_type} type alert, {event_name} \n'+
                    'Waiting 1 day to run FRA', 'blanket_blob')
        print(' - slack message sent \n')
    except Exception as e:
        print(e)
        print('Cannot post to slack (testing?)')

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
    
    #checking for events on the same day: looks for existing output files from previous runs
    count_dir=0
    for directory in os.listdir(os.environ.get('FAST_RESPONSE_OUTPUT')):
        if event_name in directory: count_dir+=1
    if count_dir==0: suffix='A'
    elif count_dir==2: suffix='B'
    elif count_dir==4: suffix='C'
    else: 
        print("COULD NOT DETERMINE EVENT SUFFIX")
        print("check for other events on the same day and re-run with args:")
        print('--skymap={} --time={} --alert_id={}'.format(skymap, str(event_mjd), run_id+':'+event_id)) 
        return
    
    print('Running {}'.format(command))
    subprocess.call([command, '--skymap={}'.format(skymap), 
        '--time={}'.format(str(event_mjd)), 
        '--alert_id={}'.format(run_id+':'+event_id),
        '--suffix={}'.format(suffix)]
        )

    event_name=event_name+suffix
    doc = False
    if args.document:
        try:
            dir_1000 = glob(os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'),
                                          '*{}_1.0e+03_s').format(event_name))
            subprocess.call([analysis_path+'document.py', '--path', dir_1000[0]])
            dir_2d = glob(os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'),
                                          '*{}_1.7e+05_s').format(event_name))
            subprocess.call([analysis_path+'document.py', '--path', dir_2d[0]])
            doc=True
        except:
            print('Failed to document to private webpage')

    try: 
        shifters = pd.read_csv(os.path.join(analysis_path,'../slack_posters/fra_shifters.csv'), 
                               parse_dates=[0,1])
        on_shift=''
        for i in shifters.index:
            if shifters['start'][i]<datetime.utcnow()<shifters['stop'][i]:
                on_shift+='<@{}> '.format(shifters['slack_id'][i])
        link = 'https://user-web.icecube.wisc.edu/~jthwaites/FastResponse/webpage/output/'
        wp_link_1000 = '{}{}_{}_1.0e+03_s.html'.format(link, eventtime[0:10].replace('-','_'),event_name)
        wp_link_2d   = '{}{}_{}_1.7e+05_s.html'.format(link, eventtime[0:10].replace('-','_'),event_name)
        bot.send_message(f'Done running FRA for {alert_type} alert, {event_name}.\n '+ on_shift +'on shift',
                         'blanket_blob')
        if doc:
            bot.send_message("Results for 1000s: <{}|link> \nResults for 2d: <{}|link>".format(
                              wp_link_1000, wp_link_2d),
                             'blanket_blob')
        print(' - slack message sent \n')
    except Exception as e:
        print(e)
        print('No slack message sent.')

if __name__ == '__main__':
    import os, subprocess
    import healpy as hp
    import numpy as np
    import lxml.etree
    import argparse
    from astropy.time import Time
    from datetime import datetime
    from dateutil.parser import parse
    import time
    from glob import glob
    from fast_response.slack_posters.slack import slackbot
    import pandas as pd
    import pwd

    username = pwd.getpwuid(os.getuid())[0]
    #default for if to document or not: only way to check reports on realtime 
    if username == 'realtime': 
        document = True
    else: 
        document = False

    parser = argparse.ArgumentParser(description='Fast Response Analysis')
    parser.add_argument('--run_live', action='store_true', default=False,
                        help='Run on live GCNs')
    parser.add_argument('--test_cascade', default=False, action='store_true',
                        help='When testing, raise to run a cascade, else track')
    parser.add_argument('--document', action='store_true', default=document,
                        help='flag to raise to push results to internal webpage')
    args = parser.parse_args()

    #for now, testing
    #with open('../slack_posters/internal_alert_slackbot.txt') as f:
    #        channel = f.readline()
    #        webhook = f.readline()
    #        bot_name = f.readline()
    #bot = slackbot(channel, bot_name, webhook)

    if args.run_live:
        print("Listening for GCNs . . . ")
        with open('../slack_posters/internal_alert_slackbot.txt') as f:
            channel = f.readline().rstrip('\n')
            webhook = f.readline().rstrip('\n')
            bot_name = f.readline().rstrip('\n')
        bot = slackbot(channel, bot_name, webhook)
        gcn.listen(handler=process_gcn)
    else:
        try:
            import fast_response
            sample_skymap_path=os.path.dirname(fast_response.__file__) +'/sample_skymaps/'
        except Exception as e:
            #future: possibly point to FRA on /data/ana/ 
            print(e)
            print('Cannot find path to sample skymaps')
            exit()
        
        if not args.test_cascade:
            print("Running on sample track . . . ")
            payload = open(sample_skymap_path \
                + 'sample_astrotrack_alert_2021.xml', 'rb').read()
            root = lxml.etree.fromstring(payload)
            process_gcn(payload, root)
        else:
            print("Running on sample cascade . . . ")
            payload = open(sample_skymap_path \
                + 'sample_cascade.txt', 'rb').read()
            root = lxml.etree.fromstring(payload)
            process_gcn(payload, root)

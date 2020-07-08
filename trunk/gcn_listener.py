''' Script to automatically receive GCN notices for IceCube
    alert events and run followup accordingly

    Author: Alex Pizzuto
    Date:   July 2020
'''

import gcn
import healpy as hp
import numpy as np
import os, subprocess
import lxml.etree

@gcn.handlers.include_notice_types(
    gcn.notice_types.ICECUBE_ASTROTRACK_GOLD,
    gcn.notice_types.ICECUBE_ASTROTRACK_BRONZE)

def process_gcn(payload, root):
    analysis_path = os.environ.get('FAST_RESPONSE_SCRIPTS')
    if analysis_path is None:
        print('###########################################################################')
        print('CANNOT FIND ENVIRONMENT VARIABLE POINTING TO REALTIME FAST RESPONSE PACKAGE\n')
        print('put \'export /path/to/realtime_gw/release\' in your bashrc')
        print('###########################################################################')
        exit()

    # Read all of the VOEvent parameters from the "What" section.
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}


    stream = params['Stream']
    event_id = params['event_id']
    run_id = params['run_id']
    pos2d = root.find('.//{*}Position2D')
    ra = float(pos2d.find('.//{*}C1').text)
    dec = float(pos2d.find('.//{*}C2').text)
    sigma = float(params['src_error_90']) / 2.14
    eventtime = root.find('.//ISOTime').text
    event_mjd = Time(eventtime, format='iso').mjd

    if 'casc' in stream.lower():
        skymap = params['skymap_fits']
        #BLA subprocess run the cascade one
    else:
        pass
    #command = analysis_path + '/run_gw_followup.py'
    #args = ['--skymap', '--trigger','--name','--role']
    #subprocess.call([command,args[0],skymap,args[1],trigger,args[2],name,args[3],root.attrib['role']])


print('Listening for GCNs...')
#gcn.listen(handler=process_gcn)

### FOR OFFLINE TESTING
payload = open('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/alert_event_followup/sample_astrotack_alert.xml', 'rb').read()
root = lxml.etree.fromstring(payload)
process_gcn(payload, root)

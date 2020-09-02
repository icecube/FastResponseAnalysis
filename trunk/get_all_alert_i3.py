import subprocess
import os
import time
from glob import glob

import healpy as hp
import numpy as np
from astropy.time import Time

skymap_files = glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run*.fits.gz')
base_output = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/alert_event_i3_files/'

for alert_ind in list(range(len(skymap_files)))[:3]:
    time.sleep(5)
    skymap_fits, skymap_header = hp.read_map(skymap_files[alert_ind], h=True, verbose=False)
    skymap_header = {name: val for name, val in skymap_header}
    run_id, ev_id = skymap_header['RUNID'], skymap_header['EVENTID']
    ev_mjd = skymap_header['EVENTMJD']
    ev_iso = skymap_header['START']
    start_iso = Time(ev_mjd - 0.00003, format='mjd').iso
    stop_iso = Time(ev_mjd + 0.00003, format='mjd').iso
    print(ev_iso)
    signalness = skymap_header['SIGNAL']
    ev_en = skymap_header['ENERGY']
    ev_ra, ev_dec = np.radians(skymap_header['RA']), np.radians(skymap_header['DEC'])
    ev_stream = skymap_header['I3TYPE']
    if 'hese' in ev_stream.lower():
        db_stream = 'HESE'
    elif 'gfu' in ev_stream.lower():
        db_stream = 'neutrino'
    elif 'ehe' in ev_stream.lower():
        db_stream = 'EHE'
    env = dict(os.environ)
    subprocess.call(['python','/data/user/apizzuto/fast_response_code/get_event_view/scripts/get_GCDQP.py',
                        '--start={}'.format(start_iso), '--stop={}'.format(stop_iso),
                        '--topic={}'.format(db_stream), 
                        '{}index_{}_run_{}_event_{}.i3'.format(base_output, alert_ind, run_id, ev_id)],
                        env = env
                       )

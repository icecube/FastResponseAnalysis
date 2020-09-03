from glob import glob
import healpy as hp

from icecube import dataclasses, dataio, icetray
from I3Tray import *

skymap_files = glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run*.fits.gz')

alert_outfile = dataio.I3File('/data/user/apizzuto/fast_response_skylab/alert_event_followup/alert_event_i3_files/concatenated_alert_evs.i3', 'w')

for alert_ind in list(range(len(skymap_files)))[:10]:
    skymap_fits, skymap_header = hp.read_map(skymap_files[alert_ind], h=True, verbose=False)
    skymap_header = {name: val for name, val in skymap_header}
    run_id, ev_id = skymap_header['RUNID'], skymap_header['EVENTID']
    ev_mjd = skymap_header['EVENTMJD']
    i3_paths = glob('/data/ana/realtime/alert_catalog_v2/i3_files/Run{:08d}_*'.format(run_id))
    if len(i3_paths) == 1:
        i3_path = i3_paths[0]
    elif ev_id == 63430929:
        i3_path = '/data/ana/realtime/alert_catalog_v2/i3_files/Run00131653_scanned1024.i3.zst'
    elif run_id == 118973:
        if ev_id == 22324184:
            i3_path = '/data/ana/realtime/alert_catalog_v2/i3_files/Run00118973_scanned1024.i3.zst'
        elif ev_id == 25391094:
            i3_path = '/data/ana/realtime/alert_catalog_v2/i3_files/Run00118973_event2_scanned1024.i3.zst'
    else:
        i3_path = glob('/data/ana/realtime/alert_catalog_v2/i3_files/Run{:08d}_event{}_*'.format(run_id, ev_id))[0]
    fi = dataio.I3File(i3_path, 'r')
    g_frame = fi.pop_frame()
    c_frame = fi.pop_frame()
    d_frame = fi.pop_frame()
    q_frame = fi.pop_frame()
    p_frame = fi.pop_frame()
    if alert_ind == 0:
        alert_outfile.push(g_frame)
        alert_outfile.push(c_frame)
        alert_outfile.push(d_frame)
    alert_outfile.push(q_frame)
    alert_outfile.push(p_frame)
    
alert_outfile.close()

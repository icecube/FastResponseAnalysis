import healpy as hp                                  
from glob import glob
import sys
sys.path.append('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/')
from FastResponseAnalysis import FastResponseAnalysis
from astropy.time import Time
TSs = []; TSDs = []
ps = []

skymap_files = glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run*.fits.gz')

for index in range(len(skymap_files)):
    print("Index: {}".format(index))
    #skymap_files = glob('/data/ana/realtime/alert_catalog_v2/2yr_prelim/fits_files/Run13*.fits.gz')
    skymap_fits, skymap_header = hp.read_map(skymap_files[index], h=True, verbose=False)
    skymap_header = {name: val for name, val in skymap_header}
    ev_mjd = skymap_header['EVENTMJD']
    if "2018-05-28" in skymap_header['START']:
        print("HUGE MAP SKIPPING")
        continue
    ev_run, ev_id = skymap_header['RUNID'], skymap_header['EVENTID']
    source = {"Skipped Events": [(ev_run, ev_id)]}
    deltaT = 2.
    event_mjd = ev_mjd
    start_mjd = event_mjd - (deltaT / 2.)
    stop_mjd = event_mjd + (deltaT / 2.)
    
    start = Time(start_mjd, format='mjd').iso
    stop = Time(stop_mjd, format='mjd').iso
    
    
    f = FastResponseAnalysis(skymap_files[index], start, stop, save=False, alert_event=True,smear=True,
                            **source)
    inj = f.initialize_injector(gamma=2.5)
    print(Time(f.start, format='mjd').iso)
    f.unblind_TS()
    TSs.append(f.ts)
    f.calc_pvalue(ntrials=100)
    ps.append(f.p)
    TSDs.append(f.tsd)
    if index % 5 == 0:
        print TSs, ps

res = {'tsd': TSDs, 'TS': TSs, 'p': ps}

import pickle
with open('/data/user/apizzuto/fast_response_skylab/dump/test_scans_subset_smear.pkl', 'w') as fi:
    pickle.dump(res, fi)

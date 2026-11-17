#!/usr/bin/env python

r"""
Run trials for background, with a given prior. 
Record TS, best-fit location, and fitted events
around best fit location

"""
import numpy as np
import healpy as hp
import os, argparse
from astropy.time import Time
import pickle

from fast_response.AlertFollowup import AlertFollowup

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--deltaT', type=float, default=None,
                    help='Time Window in seconds')
parser.add_argument('--ntrials', type=int, default = 10000,
                        help='Trials')
parser.add_argument("--alert_mjd", default=60298.0, type=float,
                        help="mjd of alert event (default Dec 20, 2023)")
parser.add_argument("--skymap", default=None, type=str,
                    help='skymap link')
parser.add_argument('--seed', default=1, type=int, 
                        help='Unique seed for running on the cluster')
parser.add_argument('--outdir',type=str, default=None,
                        help='Output directory to save npz (default = FAST_RESPONSE_OUTPUT env variable or cwd)')
args = parser.parse_args()

if args.outdir == None:
    try: 
        outdir = os.environ.get('FAST_RESPONSE_OUTPUT')
    except:
        outdir = './'
else:
    outdir = args.outdir
if not os.path.exists(outdir+'/trials/'):
    os.mkdir(outdir+'/trials/')
outdir=outdir+'/trials/'

start_mjd = args.alert_mjd - (args.deltaT / 86400. /2.)
stop_mjd = start_mjd + (args.deltaT / 86400.)
start_iso = Time(start_mjd, format='mjd').iso
stop_iso = Time(stop_mjd, format='mjd').iso
deltaT = args.deltaT / 86400.

f = AlertFollowup('Precompute_trials_test', args.skymap,
                start_iso, stop_iso, save=False)  
# f.llh.nbackground=args.bkg*args.deltaT/1000.
inj = f.initialize_injector() #just put this here to initialize f.spatial_prior
# print(f.llh.nbackground)

ntrials = args.ntrials
stop = ntrials * (args.seed+1)
start = stop-ntrials
seed_counter = start

npix = hp.nside2npix(f.nside)

ra, dec, TS_list =[], [], []
ns, seed =[], []

for jj in range(ntrials):
    # seed_counter += 1
    seed_counter = np.random.randint(0, 1000000)
    val = f.llh.scan(0.0, 0.0, scramble=True, seed = seed_counter,
            spatial_prior = f.spatial_prior, 
            time_mask = [deltaT / 2., (start_mjd + stop_mjd) / 2.],
            pixel_scan = [f.nside, f._pixel_scan_nsigma], inject = None)

    try:
        maxLoc = np.argmax(val['TS_spatial_prior_0'])  #pick out max of all likelihood ratios at diff pixels
    except ValueError:
        continue
    
    TS_list.append(val['TS_spatial_prior_0'].max())
    ns.append(val['nsignal'][maxLoc])
    ra.append(val['ra'][maxLoc])
    dec.append(val['dec'][maxLoc])
    # gamma.append(val['gamma'][maxLoc]) #fixed gamma in llh
    seed.append(seed_counter)

print("DONE")

outfilename = outdir+'seed_delta_t_{:.1e}_{}.npz'.format(args.deltaT, seed_counter)
results={
    'TS_List':TS_list,
    'ns_fit':ns,
    # 'gamma_fit':gamma,
    'ra':ra,
    'dec':dec,
    'seed':seed
}
print('saving everything')
with open(outfilename, 'wb') as f:
    pickle.dump(results, f)

print("Saved to {}".format(outfilename))
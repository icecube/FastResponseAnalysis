#!/usr/bin/env python

r"""
Run trials for background only, all-sky scans. 
No spatial prior is included; instead we record the 
TS values in each pixel and maps are applied in realtime

Modified from scripts by A. Pizzuto and R. Hussain
Jessie Thwaites, Oct 2023
"""


import numpy  as np
import healpy as hp
import argparse, os

from astropy.time           import Time
from fast_response          import GWFollowup
from skylab.llh_models      import EnergyLLH
from skylab.ps_llh          import PointSourceLLH
from skylab.spectral_models import PowerLaw 
from skylab.temporal_models import BoxProfile, TemporalModel
from scipy                  import sparse

######################### CONFIGURE ARGUEMENTS #############################
p = argparse.ArgumentParser(description="Calculates Sensitivity and Discovery"
                            " Potential Fluxes for Background Gravitational wave/Neutrino Coincidence study",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--ntrials", default=1000, type=int,
                help="Number of trials (default=1000")
p.add_argument("--seed", default=0, type=int,
                help="Process ID to save unique numpy array after running (Default=0)")
p.add_argument("--bkg", default=6.4, type=float,
                help="Expected background rate in ontime window (default=6.4)")
p.add_argument("--deltaT", default = 1000, type=int,
               help = "Time window to use, as an int! in sec (default 1000)")
p.add_argument("--outdir", default = None, type=str,
               help = "output directory")
args = p.parse_args()

def config_llh(self, scramble = True):
    """ Use full timescramble method for BG trials
    This method is adapted from initialize_llh in FastResponseAnalysis.py
    """

    self.get_data()

    print("Initializing Point Source LLH in Skylab")

    assert self._float_index==True,\
        'Index should be fitted - check initialization'

    if self._fix_index:
        llh_model = EnergyLLH(
            twodim_bins=[self.energy_bins, self.sinDec_bins],
            allow_empty=True,
            spectrum=PowerLaw(A=1, gamma=self.index, E0=1000.),
            ncpu=self._ncpu) 
    elif self._float_index:
        llh_model = EnergyLLH(
            twodim_bins=[self.energy_bins, self.sinDec_bins],
            allow_empty=True,
            bounds=self._index_range,
            seed = self.index,
            ncpu=self._ncpu)
        
    box = TemporalModel(
        grl=self.grl,
        poisson_llh=True,
        days=self._nb_days,
        signal=BoxProfile(self.start, self.stop))

    # should not be an extension here
    src_extension = None

    llh = PointSourceLLH(
        self.exp,                      # array with data events
        self.mc,                       # array with Monte Carlo events
        self.livetime,                 # total livetime of the data events
        ncpu=self._ncpu,               # use 10 CPUs when computing trials
        scramble=scramble,             # set to False for unblinding
        timescramble=True,             # not just RA scrambling
        llh_model=llh_model,           # likelihood model
        temporal_model=box,            # use box for temporal model
        nsource_bounds=(0., 1e3),      # bounds on fitted ns
        src_extension = src_extension, # Model symmetric extended source
        nsource=1.,                    # seed for nsignal fit
    )           

    return llh

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

gw_time = Time(60000., format='mjd') #2023-02-25
delta_t = float(args.deltaT)

if args.deltaT==1000:
    start_mjd = gw_time - (delta_t / 86400. / 2.)
    stop_mjd = gw_time + (delta_t / 86400. / 2.)
else: #2 week followup
    start_mjd = gw_time - 0.1
    stop_mjd = gw_time + 14.

start_iso = start_mjd.iso
stop_iso = stop_mjd.iso

#skymap needed for initialization - not used here
f = GWFollowup('precomputed_bg_test', '/data/user/jthwaites/o3-gw-skymaps/S191216ap.fits.gz', 
               start_iso, stop_iso)
f._allow_neg = False

f.llh = config_llh(f)
f.llh.nbackground=args.bkg

ntrials = args.ntrials
stop = ntrials * (args.seed+1)
start = stop-ntrials
seed_counter = start

npix = hp.nside2npix(f.nside)
shape = (args.ntrials, npix)
maps = sparse.lil_matrix(shape, dtype=float)
for jj in range(ntrials):
    val = f.llh.scan(0.0, 0.0, scramble=True, seed = seed_counter,
        #spatial_prior = f.spatial_prior, 
        time_mask = [delta_t / 2., (start_mjd.mjd + stop_mjd.mjd) / 2.],
        pixel_scan = [f.nside, f._pixel_scan_nsigma], inject = None)

    if val['TS'] is not None:
            dtype = [('ts',float),('pixel',float)]
            results = np.empty((val['TS'].size,), dtype=dtype)
            pixels = hp.ang2pix(f.nside, np.pi/2. - val['dec'], val['ra'])
            maps[jj, pixels] = val['TS']

    seed_counter+=1
print("DONE")
hp_sparse = maps.tocsr()

outfilename = outdir+'gw_{:.1f}_mHz_delta_t_{:.1e}_seed_{}.npz'.format(args.bkg, args.deltaT, args.seed)
sparse.save_npz(outfilename, hp_sparse)
print("Saved to {}".format(outfilename))
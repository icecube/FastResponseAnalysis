r"""
Config File which sets up likelihood object and PriorInjector object.
Option to choose between the two Anomalous Anita Events, namely the
event detected in the Askaryan channel ('Askaryan') or the 
event that is consistent with a tau candiate ('tau')
Modified version of Josh Wood's multi-year config file which can be found 
here: skylab/doc/analyses/icecube-170922A_wGFU/config.py 
"""


import os, sys, pickle
import numpy as np 
import healpy as hp
from time import time

from   glob                     import glob
from   skylab.ps_llh            import PointSourceLLH, MultiPointSourceLLH
from   skylab.ps_injector       import PriorInjector
from   skylab.temporal_models   import TemporalModel, BoxProfile
from   skylab.temporal_scrambling import TimeScrambler

from   skylab.llh_models        import ClassicLLH, EnergyLLH
from   skylab.datasets          import Datasets
from astropy.io                 import fits

from skylab.utils               import dist
from skylab.priors              import SpatialPrior

# energy units
GeV = 1
TeV = 1000*GeV

#@profile
def config(alert_ind, seed = 1, scramble = True, e_range=(0,np.inf), g_range=[1., 5.],
           gamma = 2.0, E0 = 1*TeV, remove = False, ncpu=20, nside=256,
           poisson=False, injector = True, verbose=True):
    r""" Configure point source likelihood and injector. 

    Parameters
    ----------
    alert_ind: int
    index of IceCube alert event

    seed : int
    Seed for random number generator

    Returns
    -------
    llh : PointSourceLLH
    Point source likelihood object
    inj : PriorInjector
     Point source injector object
    """
    seasons = [("GFUOnline_v001p01", "IC86, 2011-2018")]
    #skymaps_path = '/data/user/steinrob/millipede_scan_archive/fits_v3_prob_map/'
    #files = glob(skymaps_path + '*.fits')
    #skymap_fits = fits.open(files[alert_ind])[0].data

    #Turn this into a function read_alert_event()
    skmap_files = glob('/data/ana/realtime/alert_catalog_v2/2yr_prelim/fits_files/Run13*.fits.gz')
    skymap_f = fits.open(skymap_fies[alert_ind])
    skymap_fits = skymap_f[1].data
    skymap_header = skymap_f[1].header
    run_id, ev_id = skymap_header['RUNID'], skymap_header['EVENTID']
    ev_mjd = skymap_header['EVENTMJD']
    ev_iso = skymap_header['START']
    signalness = skymap_header['SIGNAL']
    ev_en = skymap_header['ENERGY']
    ev_ra, ev_dec = np.radians(skymap_header['RA']), np.radians(skymap_header['DEC'])
    ev_stream = skymap_header['I3TYPE']
    skymap_fits = np.exp(skymap_fits / 2.) #Convert from 2LLH to unnormalized probability

    if hp.pixelfunc.get_nside(skymap_fits)!=nside:
        skymap_fits = hp.pixelfunc.ud_grade(skymap_fits,nside)
    skymap_fits = skymap_fits/skymap_fits.sum()

    spatial_prior = SpatialPrior(skymap_fits, containment=0.99)

    llh = [] # store individual llh as lists to prevent pointer over-writing
    multillh = MultiPointSourceLLH(ncpu=1)

    if verbose:
        print("\n seasons:")
    for season in np.atleast_1d(seasons):
        sample = season[0]
        name = season[1]

        exp, mc, livetime = Datasets[sample].season(name, floor=np.radians(0.2))
        sinDec_bins = Datasets[sample].sinDec_bins(name)
        energy_bins = Datasets[sample].energy_bins(name)

        msg = "   - % 15s (" % season
        msg += "livetime %7.2f days, %6d events" % (livetime, exp.size)
        msg += ", mjd0 %.2f" % min(exp['time'])
        msg += ", mjd1 %.2f)" % max(exp['time'])
        if verbose:
            print(msg)

        llh_model = EnergyLLH(twodim_bins=[energy_bins, sinDec_bins],
                          allow_empty=True, bounds=g_range,
                          seed=gamma, kernel=1, ncpu=ncpu)

        llh.append(PointSourceLLH(exp, mc, livetime, mode="box",
                              scramble=scramble, llh_model=llh_model,
                              nsource_bounds=(0., 1e3), nsource=1.))

        multillh.add_sample(sample+ " : "+name, llh[-1])

        # save a little RAM by removing items copied into LLHs
        del exp, mc

        # END for (season)

    ######################### REMOVE EVENT

    if injector is False:
        return multillh, spatial_prior
    else:
        inj = PriorInjector(spatial_prior, seed=seed, gamma=gamma, E0=1 * TeV, bunchsize = 10)
        inj.fill(multillh.exp, multillh.mc, multillh.livetime)

        if verbose:
            print("\n injected spectrum:")
            print("   - %s" % str(inj.spectrum))

    return multillh, spatial_prior, inj


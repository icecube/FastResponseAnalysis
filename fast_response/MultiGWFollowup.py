from os.path import join
from .MultiFastResponseAnalysis import MultiPriorFollowup
from .GWFollowup import GWFollowup as GFUFollowup
from .FastResponseAnalysis import PriorFollowup

import numpy as np
import numpy.lib.recfunctions as rf
import healpy as hp

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# These methods are very similar to (Multi)PriorFollowup
# but GW analyses make some particular choices we want to reproduce here
class MultiGWFollowup(MultiPriorFollowup):

    # Defaults adopted from GWFollowup, different to PriorFollowup
    _pixel_scan_nsigma = 3.0
    _allow_neg = True
    _containment = None

    def find_coincident_events(self):
        r"""
        Find "coincident events" for a skymap
        based analysis. These are ALL ontime events,
        with a bool to indicate if they are in the 90% contour

        This works with a MultiPointSourceLLH.
        """
        if self.ts_scan is None:
            raise ValueError("Need to unblind TS before finding events")
        exp_theta = 0.5*np.pi - self.llh_exp['dec']
        exp_phi   = self.llh_exp['ra']
        exp_pix   = hp.ang2pix(self.nside, exp_theta, exp_phi)

        t_mask=(self.llh_exp['time'] <= self.stop) & (self.llh_exp['time'] >= self.start)
        events = self.llh_exp[t_mask]
        ontime_pix = hp.ang2pix(self.nside, 0.5*np.pi - events['dec'], events['ra'])
        logger.debug(f"nside={self.nside} on-time events ipix={ontime_pix}")
        overlap    = np.isin(ontime_pix, self.ipix_90)

        events = rf.append_fields(
            events, names=['in_contour', 'ts', 'ns', 'gamma', 'B'],
            data=np.empty((5, events['ra'].size)),
            usemask=False)

        for i in range(events['ra'].size):
            events['in_contour'][i]=overlap[i]
            enum_i = events['enum'][i]
            events['B'][i] = self.llh._samples[enum_i].llh_model.background(events[i])

        val_pix = self.scanned_pixels
        for i in range(events['ra'].size):
            idx, = np.where(val_pix == exp_pix[t_mask][i])
            # scan was restricted to containment fraction of spatial prior
            # if a given each event centroid overlaps, look up the scan values
            if idx.size > 0:
                events['ts'][i] = self.ts_scan['TS_spatial_prior_0'][idx[0]]
                events['ns'][i] = self.ts_scan['nsignal'][idx[0]]
                events['gamma'][i] = self.ts_scan['gamma'][idx[0]]

        self.events_rec_array = events
        self.coincident_events = [dict(zip(events.dtype.names, x)) for x  in events]
        self.save_items['coincident_events'] = self.coincident_events

    def per_event_pvalue(self):
        """
        Calculate per-event p-values. There are a few cases here: 

        - overall p < 0.1: 
            Redoes the all-sky scan, using per_event_scan, with only that single event.
            This is the same as asking the question: 
            If that single event is the only one on the sky, with this given skymap,
            what TS/p-value would we get for that event?

        - 1.0 > overall p > 0.1:
            Calculates the p-value at the reconstructed event direction. 
            Takes the TS at that location, and calculates the p-value at that location. 
            Does not re-run the scan, to save time in realtime

        - p=1.0:
            Does not get p-values for the events (all are set to None)

        """
        self.events_rec_array = rf.append_fields(
            self.events_rec_array,
            names=['pvalue'],
            data=np.empty((1, self.events_rec_array['ra'].size)),
            usemask=False
        )

        if self.tsd is None: # can not determine p-values
            for i in range(self.events_rec_array.size):
                self.events_rec_array['pvalue'][i] = None
        elif self.events_rec_array.size > 0: # always scan, most general
            for i in range(self.events_rec_array.size):
                ts, p = self.per_event_scan(self.events_rec_array[i])
                self.events_rec_array['pvalue'][i] = p
        else:
            pass # nothing to do
        
        self.coincident_events = [dict(zip(self.events_rec_array.dtype.names, x)) for x  in self.events_rec_array]
        self.save_items['coincident_events'] = self.coincident_events


class OnlyDNNOnlineFollowup(PriorFollowup):
    _base_dir = "/data/user/chraab/fast_response/multisample/extended_archival_floating"
    _sens_dir = join(_base_dir, "precomputed_sensitivity")
    _bg_dir = join(_base_dir, "precomputed_trials")
    _bg_format = '_'.join([
    'precomputed_trials_delta_t_{delta_t:.2e}',
    'nside_{nside}',
    'index_{index}',
    '{lookup}',
    'seed_*npz',
    ])
    _dataset = "DNNCascadesOnline_v001p01"
    _season_names = [f"IC86, 20{y:02d}" for y in range(24, 25+1)]
    _floor = np.radians(1.5) # can change this!
    _jitter = 3. # common default for DNN analyses
    _background_days = 60.
    _nb_days = 60.

class OnlyIceManFollowup(PriorFollowup):
    _dataset = "DNNCascadesIceMan_v001p00" 
    # limited to one season ON PURPOSE for better comparison with DNNOnline
    _season_names = [f"IC86, 20{y:02d}" for y in range(18, 18+1)]
    _jitter = 3. # common default for DNN analyses


class DNNOnlineFollowup(MultiGWFollowup):
    _base_dir = "/data/user/chraab/fast_response/multisample/extended_archival_floating"
    _sens_dir = join(_base_dir, "precomputed_sensitivity/")
    _bg_dir = join(_base_dir, "precomputed_trials/")
    _bg_format = '_'.join([
    'precomputed_trials_delta_t_{delta_t:.2e}',
    'nside_{nside}',
    'index_{index}',
    '{lookup}',
    'seed_*npz',
    ])
    _followups = [OnlyDNNOnlineFollowup]
    _fix_index = False
    _float_index = not _fix_index
    _index = 2.0
    _nside = 128

class IceManFollowup(MultiGWFollowup):
    _base_dir = "/data/user/chraab/fast_response/multisample/extended_archival_floating"
    _sens_dir = join(_base_dir, "precomputed_sensitivity/")
    _bg_dir = join(_base_dir, "precomputed_trials/")
    _bg_format = '_'.join([
    'precomputed_trials_delta_t_{delta_t:.2e}',
    'nside_{nside}',
    'index_{index}',
    '{lookup}',
    'seed_*npz',
    ])
    _followups = [OnlyIceManFollowup]
    _fix_index = False 
    _float_index = not _fix_index
    _index = 2.0
    _nside = 128

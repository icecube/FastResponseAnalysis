from os.path import join
from .MultiFastResponseAnalysis import MultiPriorFollowup
from .GWFollowup import GWFollowup as GFUFollowup
from .FastResponseAnalysis import PriorFollowup

import numpy as np

class OnlyDNNOnlineFollowup(PriorFollowup):
    _base_dir = "/data/user/chraab/fast_response/multisample/extended_archival_floating"
    _sens_dir = join(_base_dir, "precomputed_sensitivity")
    _bg_dir = join(_base_dir, "precomputed_trials")
    _bg_format = '_'.join([
    'precomputed_trials_delta_t_{delta_t:.2e}',
    'nside_{nside}',
    'index_{index}',
    '{lookup}',
    'low_stats.npz',
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


class DNNOnlineFollowup(MultiPriorFollowup):
    _base_dir = "/data/user/chraab/fast_response/multisample/extended_archival_floating"
    _sens_dir = join(_base_dir, "precomputed_sensitivity/")
    _bg_dir = join(_base_dir, "precomputed_trials/")
    _bg_format = '_'.join([
    'precomputed_trials_delta_t_{delta_t:.2e}',
    'nside_{nside}',
    'index_{index}',
    '{lookup}',
    'low_stats.npz',
    ])
    _followups = [OnlyDNNOnlineFollowup]
    _fix_index = False
    _float_index = not _fix_index
    _index = 2.0
    _nside = 128

class IceManFollowup(MultiPriorFollowup):
    _base_dir = "/data/user/chraab/fast_response/multisample/extended_archival_floating"
    _sens_dir = join(_base_dir, "precomputed_sensitivity/")
    _bg_dir = join(_base_dir, "precomputed_trials/")
    _bg_format = '_'.join([
    'precomputed_trials_delta_t_{delta_t:.2e}',
    'nside_{nside}',
    'index_{index}',
    '{lookup}',
    'low_stats.npz',
    ])
    _followups = [OnlyIceManFollowup]
    _fix_index = False 
    _float_index = not _fix_index
    _index = 2.0
    _nside = 128
    
class MultiGWFollowup(MultiPriorFollowup):
    _followups = [GFUFollowup, OnlyDNNOnlineFollowup]
    _fix_index = False
    _float_index = not _fix_index
    _index = 2.0
    
from os.path import join
from .MultiFastResponseAnalysis import MultiPriorFollowup
from .GWFollowup import GWFollowup as GFUFollowup
from .FastResponseAnalysis import PriorFollowup

import numpy as np

class OnlyDNNOnlineFollowup(PriorFollowup):
    _base_dir = "/data/user/chraab/fast_response/multisample/variable_jitter_new_environment"
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
    _season_names = ["livestream"]
    _floor = np.radians(1.5) # can change this!
    _jitter = 3. # common default for DNN analyses
    _background_days = 100.
    # Need to add analysis cuts? That's something SKATE analysers would know.

class DNNOnlineFollowup(MultiPriorFollowup):
    _base_dir = "/data/user/chraab/fast_response/multisample/variable_jitter_new_environment"
    _sens_dir = join(_base_dir, "precomputed_sensitivity")
    _bg_dir = join(_base_dir, "precomputed_trials")
    _bg_format = '_'.join([
    'precomputed_trials_delta_t_{delta_t:.2e}',
    'nside_{nside}',
    'index_{index}',
    '{lookup}',
    'low_stats.npz',
    ])
    _followups = [OnlyDNNOnlineFollowup]
    _fix_index = True # TODO check if this is what they do
    _float_index = not _fix_index
    _index = 2.0
    _nside = 64
    
class MultiGWFollowup(MultiPriorFollowup):
    _followups = [GFUFollowup, OnlyDNNOnlineFollowup]
    _fix_index = True # TODO check if this is what they do
    _float_index = not _fix_index
    _index = 2.0
    
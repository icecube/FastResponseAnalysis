
from abc import abstractmethod
import os, sys, time, subprocess
import pickle, dateutil.parser, logging, warnings
from argparse import Namespace
from copy import deepcopy

import h5py
import healpy                 as hp
import numpy                  as np
import seaborn                as sns
import matplotlib             as mpl
import matplotlib.pyplot      as plt
import numpy.lib.recfunctions as rf
from astropy.time           import Time
from scipy.special          import erfinv
from matplotlib.lines       import Line2D

from skylab.datasets        import Datasets
from skylab.llh_models      import EnergyLLH
from skylab.priors          import SpatialPrior
from skylab.ps_injector     import PointSourceInjector
from skylab.ps_llh          import PointSourceLLH, MultiPointSourceLLH
from skylab.ps_injector     import PriorInjector
from skylab.spectral_models import PowerLaw 
from skylab.temporal_models import BoxProfile, TemporalModel
# import meander

from . import web_utils
from . import sensitivity_utils
from . import plotting_utils
from .reports import FastResponseReport
from .FastResponseAnalysis import FastResponseAnalysis, PriorFollowup, PointSourceFollowup
from .FastResponseAnalysis_inherit import MultiPriorFollowup
from .GWFollowup import GWFollowup
from .AlertFollowup import AlertFollowup, CascadeFollowup, TrackFollowup

mpl.use('agg')
current_palette = sns.color_palette('colorblind', 10)
logging.getLogger().setLevel(logging.ERROR)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

# the real work comes with GWFollowup
# other than ExternalFollowup, it doesn't just add attributes
# but specify many methods, some which have to be adapted
# => the inheritance gets even deeper
# However Track and CascadeFollowup are relatively close


class MultiAlertFollowup(AlertFollowup, MultiPriorFollowup):

    # attributes that are identical for each analysis produced by one followup class
    # will be broadcast to the constituent analyses
    # this is clunky, but more explicit than some dir(...) shenanigans
    static_attributes = ['_verbose',
                         '_index',
                         '_float_index',
                         '_fix_index',
                         '_index_range',
                         '_llh_seed',
                         '_nb_days',
                         '_ncpu',
                         '_pixel_scan_nsigma',
                         '_allow_neg',
                         '_containment',
                         '_nside',
                         '_smear',
                         ]
    _base_dir = '/data/user/chraab/fast_response/multisample_alert_followup/'
    _bg_trial_dir = os.path.join(_base_dir, 'alert_precomputed_trials')
    _sens_dir = os.path.join(_base_dir, 'reference_sensitivity_curves')

    def run_background_trials(self, ntrials = 1000):
        # How to get back to base class method?
        raise NotImplementedError('')
    
    # def run_background_trials(self, ntrials = 1000):
        # TODO  this is tricky - don't implement caching for now.
        # as it deals with Multi-LLH trials (calculated where?)
        # but uses self.llh.nbackground (which is defined for PSLLH only)
        # to retrieve BG trials for the current rate - in multi
        # and then it also requires a strategy for background rate in multiple samples.
        # raise NotImplementedError('''Not decided how to handle background fluctuations
        #                           in multiple samples, which may be partially correlated,
        #                           but energy x selection adds another axis to seasonality.''')
        #duration_seconds = self.duration * 86400.
        # pre_ts_array = sparse.load_npz(os.path.join(
        #     bg_trial_dir,
        #     f'precomputed_trials_delta_t_{duration_seconds:.2e}_low_stats.npz',
        #)
        # ts_norm = np.log(np.amax(self.skymap))
        # ts_prior = pre_ts_array.copy()
        # ts_prior.data += 2.*(np.log(self.skymap[pre_ts_array.indices]) - ts_norm)
        # ts_prior.data[~np.isfinite(ts_prior.data)] = 0.
        # ts_prior.data[ts_prior.data < 0] = 0.
        # tsd = ts_prior.max(axis=1).A
        # tsd = np.array(tsd)
        # self.tsd = tsd
        # return tsd

class MultiTrackFollowup(MultiAlertFollowup, TrackFollowup):
    # this inheritance structure should work as long as CascadeFollowup does not override methods
    # but getting very spaghetti already and some refactoring might be in order
    _smear = True

    
class MultiCascadeFollowup(MultiAlertFollowup, CascadeFollowup):
    _smear = False

    
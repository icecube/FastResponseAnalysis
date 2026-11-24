from .MultiFastResponseAnalysis import MultiPointSourceFollowup
from .FastResponseAnalysis import PointSourceFollowup
from .ExternalFollowup import ExternalFollowup as GFUFollowup

import numpy as np

# Consistent with existing structure:
# base class: methods
# specific class: data sample, class attributes, "analysis definition"
# instance: specific source/follow-up
class GrecoFollowup(PointSourceFollowup):
    """
    An externally-triggered point source followup
    compatible to the GFU-based ExternalFollowup
    but using the GrecoOnline dataset and respective settings,
    currently only archival data and used for sensitivity studies.
    """
    _dataset = 'GrecoOnline_v002p09'
    _season_names = [f"IC86, 20{y:02d}" for y in range(12, 19 + 1)]
    _floor = np.radians(0.2) # can change this!

class DNNFollowup(PointSourceFollowup):
    """
    An externally-triggered point source followup
    compatible to the GFU-based ExternalFollowup
    but using the archival DNNCascades dataset and respective settings,
    used for sensitivity studies.
    """
    _dataset = "DNNCascades_v001p01"
    _season_names = [f"IC86, 20{y:02d}" for y in range(18, 19+1)]
    _floor = np.radians(1.5) # can change this!
    _jitter = 3. # common default for DNN analyses
    
    
class DNNOnlineFollowup(PointSourceFollowup):
    """
    An externally-triggered point source followup
    compatible to the GFU-based ExternalFollowup
    but using the DNNCascadesOnline dataset and respective settings,
    currently not yet approved for analysis.
    """
    _dataset = "DNNCascadesOnline_v001p01"
    _season_names = [f"IC86, 20{y:02d}" for y in range(24, 25+1)]
    _floor = np.radians(1.5) # can change this!
    _jitter = 3. # common default for DNN analyses
    # Background loading / calculation windows:
    # about x10 vs. GFU to get good rate estimate for 1/30th of the rate
    _background_days = 60.
    _nb_days = 60.
    

class DNNIceManFollowup(PointSourceFollowup):
    """
    An externally-triggered point source followup
    compatible to the GFU-based ExternalFollowup
    but using the archival DNNCascades-IceMan dataset,
    currently not yet approved for analysis.
    """
    _dataset = "DNNCascadesIceMan_v001p00" 
    # limited to one season ON PURPOSE for better comparison with DNNOnline
    _season_names = [f"IC86, 20{y:02d}" for y in range(18, 18+1)]
    _jitter = 3. # common default for DNN analyses

class MultiFollowup(MultiPointSourceFollowup):
    """
    Class for external point-source or extended source followup.
    By default, uses floating index of 2.5 in the LLH. Based on 
    the PointSourceFollowup class adapted to accept multiple samples
    via the configuration of constituent follow-up analyses.
    The floating spectral index is broadcast to the constituent datasets/likelihoods.
    """
    _followups = [GFUFollowup, GrecoFollowup]
    _fix_index = False # => this 
    _float_index = not _fix_index
    _index = 2.5
    
    
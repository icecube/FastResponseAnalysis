from .MultiFastResponseAnalysis import MultiPointSourceFollowup
from .FastResponseAnalysis import PointSourceFollowup

import numpy as np

# Consistent with existing structure:
# base class: methods
# specific class: data sample, class attributes, "analysis definition"
# instance: specific source/follow-up
class GFUFollowup(PointSourceFollowup):
    _dataset = "GFUOnline_v001p02"
    _season_names = [f"IC86, 201{y}" for y in range(7, 10)]
    #_season_names = [f"IC86, 201{y}" for y in range(1, 10)]
    _floor = np.radians(0.2)

class GrecoFollowup(PointSourceFollowup):
    _dataset = 'GrecoOnline_v002pFactor'
    _season_names = [f"IC86, 20{y:02d}" for y in range(18, 20+1)]
    #_season_names = [f"IC86, 20{y:02d}" for y in range(12, 22+1)]
    _floor = np.radians(0.2) # can change this!
    # extended/GRB-style LLH?

class DNNFollowup(PointSourceFollowup):
    _dataset = "DNNCascades_v001p01" # TODO switch to "online" version
    _season_names = [f"IC86, 20{y:02d}" for y in range(18, 19+1)]
    _floor = np.radians(1.5) # can change this!
    _jitter = 3. # common default for DNN analyses
    
    
#or any other number of definitions

class DNNOnlineFollowup(PointSourceFollowup):
    _dataset = "DNNCascadesOnline_v001p01"
    _season_names = [f"IC86, 20{y:02d}" for y in range(24, 25+1)]
    _floor = np.radians(1.5) # can change this!
    _jitter = 3. # common default for DNN analyses
    _background_days = 60.
    _nb_days = 60.
    

class DNNIceManFollowup(PointSourceFollowup):
    _dataset = "DNNCascadesIceMan_v001p00" 
    # limited to one season ON PURPOSE for better comparison with DNNOnline
    _season_names = [f"IC86, 20{y:02d}" for y in range(18, 18+1)]
    _jitter = 3. # common default for DNN analyses

class MultiFollowup(MultiPointSourceFollowup):
    '''
    Class for external point-source or extended source followup.
    By default, uses floating index of 2.5 in the LLH. Based on 
    the PointSourceFollowup class adapted to accept multiple samples
    via the configuration of constituent follow-up analyses.

    '''
    _followups = [GFUFollowup, GrecoFollowup] # more consistent
    _fix_index = False
    _float_index = not _fix_index
    _index = 2.5
    
    
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
    
#or any other number of definitions

class MultiFollowup(MultiPointSourceFollowup):
    '''
    Class for external point-source or extended source followup.
    By default, uses floating index of 2.5 in the LLH. Based on 
    the PointSourceFollowup class adapted to accept multiple samples
    via the configuration of constituent follow-up analyses.

    '''
    _followups = [GFUFollowup, GrecoFollowup] # more consistent
    _fix_index = True
    _float_index = not _fix_index
    _index = 2.5
    
    
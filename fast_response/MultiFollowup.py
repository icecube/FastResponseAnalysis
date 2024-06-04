from .FastResponseAnalysis_ifelse import PointSourceFollowup

class ExternalFollowup(PointSourceFollowup):
    '''
    Class for external point-source or extended source followup.
    By default, uses a fixed index of 2.0 in LLH. Based on 
    the PointSourceFollowup class

    '''
    _dataset = ["GFUOnline_v001p02",
                "GrecoOnline_v002p10",
               ]
    
    _fix_index = True
    _float_index = not _fix_index
    _index = 2.5
    _scramble = True
from .FastResponseAnalysis import PointSourceFollowup, PriorFollowup

class ExternalFollowup(PointSourceFollowup):
    '''
    Class for external point-source or extended source followup.
    By default, uses a fixed index of 2.0 in LLH. Based on 
    the PointSourceFollowup class

    '''
    _dataset = "GFUOnline_v001p02"
    _fix_index = True
    _float_index = not _fix_index
    _index = 2.0

class ExternalSkymapFollowup(PriorFollowup):
    '''
    Class for external skymap followup.
    By default, uses a fixed index of 2.0 in LLH.
    Based on the PriorFollowup class.

    NOTE: There are many functions that are not yet implemented here. 
    The point source followup is strongly preferred for this reason.
    '''
    _dataset = "GFUOnline_v001p02"
    _fix_index = True
    _float_index = not _fix_index
    _index = 2.0

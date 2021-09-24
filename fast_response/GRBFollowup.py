from .FastResponseAnalysis import PriorFollowup

class GRBFollowup(PriorFollowup):
    _dataset = 'GFUOnline_v001_p03'
    _fix_index = True
    _index = 2.0

    pass

    # index fixed to 2
    # Dataset: with zenith_smoothing
    # Different precomputed trials
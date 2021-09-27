from FastResponseAnalysis import PointSourceFollowup

class ExternalFollowup(PointSourceFollowup):
    _dataset = "GFUOnline_v001p02"
    _fix_index = True
    _float_index = not _fix_index
    _index = 2.0
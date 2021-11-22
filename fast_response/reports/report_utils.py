from astropy.time import Time, TimeDelta
import numpy as np

def interptime(tstr):
    if tstr is None:
        t= None
    elif tstr.endswith('sec'):
        t = TimeDelta(float(tstr[:-3]),format='sec')
    elif tstr.endswith('s'):
        t = TimeDelta(float(tstr[:-1]),format='sec')
    elif tstr.endswith('min'):
        t = TimeDelta(float(tstr[:-3])*60,format='sec')
    elif tstr.endswith('hr'):
        t = TimeDelta(float(tstr[:-2])*3600,format='sec')
    elif tstr.startswith('mjd'):
        t = Time(float(tstr[3:]),format='mjd')
    elif tstr.endswith('d'):
        t = TimeDelta(float(tstr[:-1]),format='jd')
    else:
        t = Time(tstr,format='iso')
    return t
 
def get_times(trig,sta,sto):
    trigger = interptime(trig)
    if type(trigger)!=Time:
        trigger=None
    start = interptime(sta)

    if type(trigger)==Time and type(start)==Time:
        pass
    elif type(trigger)==Time and start is None:
        start = trigger
    elif type(start)==Time and trigger is None:
        trigger = start
    elif type(start)==TimeDelta and trigger is not None:
        start = trigger+start
    else:
        raise Exception("either --trigger or --start must be a absolute time")

    stop=interptime(sto)
    if type(stop)!=Time:
        stop = trigger+stop

    return trigger,start,stop

def valid_location(ra, dec):
    r'''For a point source, check if location
    is valid
    Parameters:
    -----------
    ra: float
        Right ascension in radians
    dec: float
        Declination in radians
    Returns:
    --------
    True or False: Bool
        Whether or not the location is valid
    '''
    if (ra is None) or (dec is None):
        return False
    else:
        if (ra < 0.) or (ra > 2*np.pi):
            return False
        elif (dec < -np.pi / 2.) or (dec > np.pi/2.):
            return False
        else:
            return True
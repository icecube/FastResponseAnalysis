import scipy as sp
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def find_nearest(array, value):
    array = np.asarray(array)
    ind = (np.abs(array - value)).argmin()
    return array[ind]

def deltaPsi(dec1, ra1, dec2, ra2):
    """
    Calculate angular distance.
    
    Args:
        dec1: Declination of first direction in radian
        ra1: Right ascension of first direction in radian
        dec2: Declination of second direction in radian
        ra2: Right ascension of second direction in radian
        
    Returns angular distance in radian
    """
    return deltaPsi2(np.sin(dec1), np.cos(dec1), np.sin(ra1), np.cos(ra1), np.sin(dec2), np.cos(dec2), np.sin(ra2), np.cos(ra2))

def deltaPsi2(sDec1, cDec1, sRa1, cRa1, sDec2, cDec2, sRa2, cRa2):
    """
    Calculate angular distance.
    
    Args:
        sDec1: sin(Declination of first direction)
        cDec1: cos(Declination of first direction)
        sRa1: sin(Right ascension of first direction)
        cRa1: cos(Right ascension of first direction)
        sDec2: sin(Declination of second direction)
        cDec2: cos(Declination of second direction)
        sRa2: sin(Right ascension of second direction)
        cRa2: cos(Right ascension of second direction)
        
    Returns angular distance in radian
    """
    tmp = cDec1*cRa1*cDec2*cRa2 + cDec1*sRa1*cDec2*sRa2 + sDec1*sDec2
    tmp[tmp>1.] = 1.
    tmp[tmp<-1.] = -1.
    return np.arccos(tmp)

def sensitivity_fit(signal_fluxes, passing, errs, fit_func, p0 = None, conf_lev = 0.9):
    r'''
    Given an array of injected fluxes (or events) and passing fraction
    relative to the unblinded TS, do a fit given a CDF
    '''
    try:
        name = fit_func.__name__
        name = name.replace("_", " ")
    except:
        name = 'fit'
    popt, pcov = curve_fit(fit_func, signal_fluxes, passing, sigma = errs, p0 = p0, maxfev=50000)
    fit_points = fit_func(signal_fluxes, *popt)
    chi2 = np.sum((fit_points - passing)**2. / errs**2.)
    dof = len(fit_points) - len(popt)
    xfit = np.linspace(np.min(signal_fluxes) - 0.5, np.max(signal_fluxes), 100)
    yfit = fit_func(xfit, *popt)
    pval = sp.stats.chi2.sf(chi2, dof)
    sens = xfit[find_nearest_idx(yfit, 0.9)]
    return {'popt': popt, 'pcov': pcov, 'chi2': chi2, 
            'dof': dof, 'xfit': xfit, 'yfit': yfit, 
            'name': name, 'pval':pval, 'ls':'--', 'sens': sens}

def erfunc(x, a, b):
    x = np.array(x)
    return 0.5 + 0.5*sp.special.erf(a*x + b)

def chi2cdf(x, df1, loc, scale):
    func = chi2.cdf(x,df1,loc,scale)
    return func

def incomplete_gamma(x, a, scale):
    x = np.array(x)
    return sp.special.gammaincc( scale*x, a)

def poissoncdf(x, mu, loc):
    func = sp.stats.poisson.cdf(x, mu, loc)
    return func

def binomial_error(p, number):
    errs = np.sqrt(p*(1.-p) / number)
    ntrig = p * number
    bound_case_pass = (ntrig + (1./3.)) / (number + (2./3.))
    bound_case_sigma = np.sqrt(bound_case_pass*(1. - bound_case_pass) / (number + 2))
    errs = np.maximum(errs, bound_case_sigma)
    return errs
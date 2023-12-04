import scipy as sp
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2

def find_nearest_idx(array, value):
    """
    Find the index where an array is nearest to a 
    given value. For the closest value, rather than index,
    see find_nearest.
    
    Parameters
    -----------
    array: list or array
        Array to check
    value: float
        Particular value to look for in an array

    Returns
    ----------
    int: 
        Index where array is closest to value

    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def find_nearest(array, value):
    """
    Find the value of an array closest to a 
    given value. For the index, rather than value, see
    find_nearest_idx.
    
    Parameters
    -----------
    array: list or array
        Array to check
    value: float
        Particular value to look for in an array

    Returns
    ----------
    float: 
        value of the array that is closest to the value given

    """
    array = np.asarray(array)
    ind = (np.abs(array - value)).argmin()
    return array[ind]

def deltaPsi(dec1, ra1, dec2, ra2):
    """
    Calculate angular distance between two given points 
    Values in Right Ascension/declination 
    
    Parameters
    -----------
    dec1: float
        Declination of first direction in radian
    ra1: float
        Right ascension of first direction in radian
    dec2: float
        Declination of second direction in radian
    ra2: float
        Right ascension of second direction in radian
    
    Returns
    -----------
    float: 
        angular distance in radian

    """
    return deltaPsi2(np.sin(dec1), np.cos(dec1), np.sin(ra1), np.cos(ra1), np.sin(dec2), np.cos(dec2), np.sin(ra2), np.cos(ra2))

def deltaPsi2(sDec1, cDec1, sRa1, cRa1, sDec2, cDec2, sRa2, cRa2):
    """
    Helper function to calculate angular distance. 
    Should not be called directly, rather call deltaPsi which
    calls this helper function.
    
    Parameters
    -----------
    sDec1: float
        sin(Declination of first direction)
    cDec1: float
        cos(Declination of first direction)
    sRa1: float
        sin(Right ascension of first direction)
    cRa1: float
        cos(Right ascension of first direction)
    sDec2: float
        sin(Declination of second direction)
    cDec2: float
        cos(Declination of second direction)
    sRa2: float
        sin(Right ascension of second direction)
    cRa2: float
        cos(Right ascension of second direction)
    
    Returns
    ---------
    float: 
        angular distance in radian

    """
    tmp = cDec1*cRa1*cDec2*cRa2 + cDec1*sRa1*cDec2*sRa2 + sDec1*sDec2
    tmp[tmp>1.] = 1.
    tmp[tmp<-1.] = -1.
    return np.arccos(tmp)

def sensitivity_fit(signal_fluxes, passing, errs, fit_func, p0 = None, conf_lev = 0.9):
    r"""
    Given an array of injected fluxes (or events) and passing fraction
    relative to the unblinded TS, do a fit given a CDF. 
    Built on curve_fit from scipy.optimize.

    Parameters
    -----------
    signal_fluxes: float Numpy array
        set of injected flux values (or signal events).
    passing: float Numpy array
        corresponding fraction of trials with TS > threshold TS value, for each signal flux value.
    errs: float Numpy array
        binomial error for the given passing values. See also: binomial_error
    fit_func: function
        function to use in the fit. 
        See also: Error function (erfunc), Chi squared (chi2cdf), Incomplete gamma (incomplete_gamma)
    p0: None, scalar or N-length sequence
        initial guess for parameters, passed directly to scipy.optimize.curve_fit
    conf_lev: float
        confidence level used (default 0.9 = 90% CL)
    
    Returns
    ----------
    dict:
        Fit values, with keys:
         
        * popt, pcov: returned fit values from scipy.optimize.curve_fit
        * chi2, dof, pval: chi squared, degrees of freedom, and p-value from fit
        * xfit, yfit: arrays of values, using the fitted parameters of the given function
        * name: name of fit used
        * ls: linestyle, returns --
        * sens: calculated sensitivity value, in flux or ns (whichever is passed in signal_fluxes)
    
    """
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
    sens = xfit[find_nearest_idx(yfit, conf_lev)]
    return {'popt': popt, 'pcov': pcov, 'chi2': chi2, 
            'dof': dof, 'xfit': xfit, 'yfit': yfit, 
            'name': name, 'pval':pval, 'ls':'--', 'sens': sens}

def erfunc(x, a, b):
    """Error function"""
    x = np.array(x)
    return 0.5 + 0.5*sp.special.erf(a*x + b)

def chi2cdf(x, df1, loc, scale):
    """Chi-squared CDF"""
    func = chi2.cdf(x,df1,loc,scale)
    return func

def incomplete_gamma(x, a, scale):
    """Incomplete gamma function"""
    x = np.array(x)
    return sp.special.gammaincc( scale*x, a)

def poissoncdf(x, mu, loc):
    """Poisson CDF function"""
    func = sp.stats.poisson.cdf(x, mu, loc)
    return func

def binomial_error(p, number):
    """
    Calculates binomial errors for given passing fraction values

    Parameters
    -----------
    p: float Numpy array
        passing fraction values
    number: int Numpy array
        number of signal trials run for each injected value
    
    Returns
    ----------
    float array:
        Numpy array of calculated binomial errors
    """
    errs = np.sqrt(p*(1.-p) / number)
    ntrig = p * number
    bound_case_pass = (ntrig + (1./3.)) / (number + (2./3.))
    bound_case_sigma = np.sqrt(bound_case_pass*(1. - bound_case_pass) / (number + 2))
    errs = np.maximum(errs, bound_case_sigma)
    return errs
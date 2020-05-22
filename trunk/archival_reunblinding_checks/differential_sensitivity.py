import numpy as np
import scipy as sp
from glob import glob
import os, sys, pickle, argparse
from astropy.time import Time
from astropy.time import TimeDelta
from scipy.optimize       import curve_fit
from scipy.stats          import chi2
from numpy.lib.recfunctions import append_fields
sys.path.append('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/')
from FastResponseAnalysis import FastResponseAnalysis

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--deltaT', type=float, default=None,
                    help='Time Window in seconds')
parser.add_argument('--sinDec', type=float, default=None,
                    help='sinDec of source')
args = parser.parse_args()

def calc_errs(signal_fluxes, passing, ntr=20):
    errs = np.sqrt(passing*(1. - passing) / ntr)
    ntrig = passing * ntr
    bound_case_pass = (ntrig + (1./3.)) / (ntr + (2./3.))
    bound_case_sigma = np.sqrt(bound_case_pass*(1. - bound_case_pass) / (ntr + 2))
    errs = np.maximum(errs, bound_case_sigma)
    return signal_fluxes, passing, errs

def sensitivity_fit(signal_fluxes, passing, errs, fit_func, p0 = None, conf_lev = 0.9):
    try:
        name = fit_func.__name__
        name = name.replace("_", " ")
    except:
        name = 'fit'
    popt, pcov = curve_fit(fit_func, signal_fluxes, passing, sigma = errs, p0 = p0, maxfev=10000)
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

def calc_sensitivity(signal_fluxes, passing, threshold = 0.5, in_ns = True, with_err = True, 
        trim=-1, p0=None, conf_lev=0.9):
    signal_fluxes, passing, errs = calc_errs(signal_fluxes, passing)
    fits, plist = [], []
    try:
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, chi2cdf, p0=p0, conf_lev=conf_lev))
        plist.append(fits[-1]['pval'])
    except:
        pass
    try:
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, erfunc, p0=p0, conf_lev=conf_lev))
        plist.append(fits[-1]['pval'])
    except:
        pass
    try:
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, incomplete_gamma, p0=p0, conf_lev=conf_lev))
        plist.append(fits[-1]['pval'])
    except:
        pass
    #Find best fit of the three, make it look different in plot
    plist = np.array(plist)
    best_fit_ind= np.argmax(plist)
    return fits[best_fit_ind]

def erfunc(x, a, b):
    return 0.5 + 0.5*sp.special.erf(a*x + b)

def chi2cdf(x,df1,loc,scale):
    func = chi2.cdf(x,df1,loc,scale)
    return func

def incomplete_gamma(x, a, scale):
    x = np.array(x)
    return sp.special.gammaincc( scale*x, a)

def poissoncdf(x, mu, loc):
    func = sp.stats.poisson.cdf(x, mu, loc)
    return func

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


gamma = 2.0
nsigs = np.array([1., 1.5, 2., 2.5, 3., 4., 5., 6., 8.])
deltaT = args.deltaT / 86400.
sinDec = args.sinDec
dec = np.arcsin(sinDec)
start_mjd = 57874.0
start = Time(start_mjd, format='mjd').iso
stop = (Time(start) + TimeDelta(deltaT)).iso
sensitivity, discovery = [], []
sensitivity_events, discovery_events = [], []
ntrials = 100
if sinDec < -0.1:
    e_bins = np.logspace(4., 8., 5)
else:
    e_bins = np.logspace(2., 8., 7)

f = FastResponseAnalysis("0., {}".format(dec*180. / np.pi), start, stop, save=False, alert_event=True)
print("Running Background")
bg = f.llh.do_trials(3000, src_ra = 0., src_dec = dec, injector = None) #, mean_signal=0.0)
bg_med = np.percentile(bg['TS'], 50.)

for lowE in e_bins:
    sens_passing = []
    fluxes = []
    #f = FastResponseAnalysis("0., {}".format(dec*180. / np.pi), start, stop, save=False, alert_event=True)
    inj = f.initialize_injector(gamma=gamma, e_range=(lowE, lowE*10.))
    #print("Running Background")
    #bg = f.llh.do_trials(3000, src_ra = 0., src_dec = dec, injector = inj, mean_signal=0.0)
    #bg_med = np.percentile(bg['TS'], 50.)
    print("Starting signal injection")
    for nsig in nsigs:
        result = f.llh.do_trials(ntrials, src_ra = 0., src_dec = dec, injector = inj, mean_signal=nsig, poisson=True)
        fluxes.append(inj.mu2flux(int(nsig)))
        sens_passing.append(np.count_nonzero(result['TS'] > bg_med))

    sens_passing = np.array(sens_passing) / float(ntrials)

    print(sens_passing)
    print(calc_errs(nsigs, sens_passing))
    sensitivity.append(fluxes[0] * calc_sensitivity(nsigs, sens_passing)['sens'])
    sensitivity_events.append(calc_sensitivity(nsigs, sens_passing)['sens'])

results = {'sensitivity': sensitivity,  
           'low_energies': e_bins,
            'high_energies': e_bins*10., 
          'sens_evs': sensitivity_events}

with open('/data/user/apizzuto/fast_response_skylab/dump/differential_sensitivity_deltaT_{:.2e}_sinDec_{:.2f}.pkl'.format(deltaT, sinDec), 'wb') as f:
    pickle.dump(results, f)

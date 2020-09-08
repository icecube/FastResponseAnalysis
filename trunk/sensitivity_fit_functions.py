import numpy as np
import scipy as sp
#from lmfit                import Model
from scipy.optimize       import curve_fit
from scipy.stats          import chi2
import pickle
from glob import glob
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/time_integrated_scripts/')
import steady_sensitivity_fits
import healpy as hp
 
palette = ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f']
skymap_files = glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run*.fits.gz')
l_ind = skymap_files[0].find("Run")
r_ind = skymap_files[0].find("_nside")

def erfunc(x, a, b):
    return 0.5 + 0.5*sp.special.erf(a*x + b)

def chi2cdf(x,df1,loc,scale):
    func = chi2.cdf(x,df1,loc,scale)
    return func

def fsigmoid(x, a, b):
    return 1.0 / (1.0 + np.exp(-a*(x-b)))

def incomplete_gamma(x, a, scale):
    return sp.special.gammaincc( scale*x, a)

def poissoncdf(x, mu, loc):
    func = sp.stats.poisson.cdf(x, mu, loc)
    return func

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def n_to_flux(N, index, delta_t, smear=True):
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/sensitivity/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    with open(fs[0], 'r') as f:
        signal_trials = pickle.load(f)
    fl_per_one = np.mean(np.array(signal_trials['flux']) / np.array(signal_trials['mean_ninj']))
    return fl_per_one * N

def dec_of_map(index, smear=True):
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/sensitivity/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, 1000.0))
    with open(fs[0], 'r') as f:
        signal_trials = pickle.load(f)
    ra, dec = np.median(signal_trials['ra']), np.median(signal_trials['dec'])
    return ra, dec

def background_distribution(index, delta_t, smear=True):
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/bg/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    with open(fs[0], 'r') as f:
        bg_trials = pickle.load(f)
    return bg_trials['ts_prior']

def signal_distribution(index, delta_t, ns, smear=True):
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/sensitivity/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    with open(fs[0], 'r') as f:
        signal_trials = pickle.load(f)
    ret = {}
    msk = np.array(signal_trials['mean_ninj']) == ns
    for k, v in signal_trials.iteritems():
        ret[k] = np.array(v)[msk]
    return ret

def pass_vs_inj(index, delta_t, threshold = 0.5, in_ns = True, with_err = True, trim=-1, smear=True):
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/bg/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    with open(fs[0], 'r') as f:
        bg_trials = pickle.load(f)
    bg_trials = bg_trials['ts_prior']
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/sensitivity/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    with open(fs[0], 'r') as f:
        signal_trials = pickle.load(f)
    bg_thresh = np.percentile(bg_trials, threshold * 100.)
    #print(bg_thresh)
    signal_fluxes, signal_indices = np.unique(signal_trials['mean_ninj'], return_index=True)
    signal_indices = np.append(signal_indices, len(signal_trials['ts_prior']))
    #print signal_indices, signal_fluxes
    if trim != -1 and trim < 0:
        signal_indices = signal_indices[:trim]
        signal_fluxes = signal_fluxes[:trim]
    elif trim > 0:
        signal_indices = signal_indices[:trim + 1]
        signal_fluxes = signal_fluxes[:trim]
    #print signal_trials['ts_prior'], signal_trials['mean_ninj']
    passing = np.array([np.count_nonzero(signal_trials['ts_prior'][li:ri] > bg_thresh) / float(ri - li) for li, ri in zip(signal_indices[:-1], signal_indices[1:])])
    if not with_err:
        return signal_fluxes, passing
    else:
        errs = np.array([np.sqrt(p*(1.-p) / float(ri - li)) for p, li, ri in zip(passing, signal_indices[:-1], signal_indices[1:])])
        ngen = np.array([float(ri - li) for li, ri in zip(signal_indices[:-1], signal_indices[1:])])
        ntrig = passing * ngen
        bound_case_pass = (ntrig + (1./3.)) / (ngen + (2./3.))
        bound_case_sigma = np.sqrt(bound_case_pass*(1. - bound_case_pass) / (ngen + 2))
        errs = np.maximum(errs, bound_case_sigma)
        return signal_fluxes, passing, errs
    
def sensitivity_curve(index, delta_t, threshold = 0.5, in_ns = True, with_err = True, trim=-1, ax = None, 
                    p0 = None, fontsize = 16, conf_lev = 0.9, smear=True, legend=True, text=True):
    signal_fluxes, passing, errs = pass_vs_inj(index, delta_t, threshold=threshold, in_ns=in_ns, with_err=with_err, trim=trim, smear=smear)
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
    try:
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, fsigmoid, p0=p0, conf_lev=conf_lev))
        plist.append(fits[-1]['pval'])
    except:
        pass
        #print("at least one fit failed")
    #Find best fit of the three, make it look different in plot
    plist = np.array(plist)
    if len(plist) > 0:
        best_fit_ind= np.argmax(plist)
        if fits[best_fit_ind]['chi2'] / fits[best_fit_ind]['dof'] < 5:
            fits[best_fit_ind]['ls'] = '-'
    
    if ax==None:
        fig, ax = plt.subplots()
    
    for fit_dict in fits:
        ax.plot(fit_dict['xfit'], fit_dict['yfit'], 
                 label = r'{}: $\chi^2$ = {:.2f}, d.o.f. = {}'.format(fit_dict['name'], fit_dict['chi2'], fit_dict['dof']),
                ls = fit_dict['ls'])
        if fit_dict['ls'] == '-':
            ax.axhline(conf_lev, color = palette[-1], linewidth = 0.3, linestyle = '-.')
            ax.axvline(fit_dict['sens'], color = palette[-1], linewidth = 0.3, linestyle = '-.')
            if text:
                ax.text(6, 0.5, r'Sens. = {:.2f}'.format(fit_dict['sens']))
    if fits[best_fit_ind]['chi2'] / fits[best_fit_ind]['dof'] > 5:
        inter = np.interp(conf_lev, passing, signal_fluxes)
        ax.axhline(conf_lev, color = palette[-1], linewidth = 0.3, linestyle = '-.')
        ax.axvline(inter, color = palette[-1], linewidth = 0.3, linestyle = '-.')
    ax.errorbar(signal_fluxes, passing, yerr=errs, capsize = 3, linestyle='', marker = 's', markersize = 2)
    if legend:
        ax.legend(loc=4, fontsize = fontsize)
    
def calc_sensitivity(index, delta_t, threshold = 0.5, in_ns = True, with_err = True, trim=-1, 
                    conf_lev = 0.9, p0=None, smear=True):
    signal_fluxes, passing, errs = pass_vs_inj(index, delta_t, threshold=threshold, in_ns=in_ns, with_err=with_err, trim=trim, smear=smear)
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
    try:
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, fsigmoid, p0=p0, conf_lev=conf_lev))
        plist.append(fits[-1]['pval'])
    except:
        pass
    #Find best fit of the three, make it look different in plot
    plist = np.array(plist)
    if len(plist) > 0:
        best_fit_ind= np.argmax(plist)
        if fits[best_fit_ind]['chi2'] / fits[best_fit_ind]['dof'] < 5:
            return fits[best_fit_ind]
    inter = np.interp(conf_lev, passing, signal_fluxes)
    return {'sens': inter, 'name': 'linear_interpolation'}
    
def sensitivity_fit(signal_fluxes, passing, errs, fit_func, p0 = None, conf_lev = 0.9):
    try:
        name = fit_func.__name__
        name = name.replace("_", " ")
    except:
        name = 'fit'
    popt, pcov = curve_fit(fit_func, signal_fluxes, passing, sigma = errs, p0 = p0, maxfev=10000)
    #print popt
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

def pvals_for_signal(index, delta_t, ns = 1, sigma_units = False, smear=True):
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/bg/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    with open(fs[0], 'r') as f:
        bg_trials = pickle.load(f)
    bg_trials = bg_trials['ts_prior']
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/sensitivity/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    with open(fs[0], 'r') as f:
        signal_trials = pickle.load(f)
    pvals = [100. - sp.stats.percentileofscore(bg_trials, ts, kind='strict') for ts in signal_trials['ts_prior']]
    pvals = np.array(pvals)*0.01
    pvals = np.where(pvals==0, 1e-6, pvals)
    if not sigma_units:
        return pvals
    else:
        return sp.stats.norm.ppf(1. - (pvals / 2.))

def find_all_sens(delta_t, smear=True, with_disc=True, disc_conf=0.5, 
                    disc_thresh=1.-0.0013, verbose=True):
    num_alerts = 249
    sensitivities = np.zeros(num_alerts)
    if with_disc:
        discoveries = np.zeros(num_alerts)
    for ind in range(num_alerts):
        if verbose:
            print ind, 
        try:
            sens = n_to_flux(calc_sensitivity(ind, delta_t, smear=smear)['sens'], 
                                ind, delta_t, smear=smear)
            sensitivities[ind] = sens
            if with_disc:
                disc = n_to_flux(calc_sensitivity(ind, delta_t, threshold=disc_thresh, 
                                    conf_lev=disc_conf, smear=smear)['sens'], ind, delta_t, smear=smear)
                discoveries[ind] = disc
            if sens*delta_t*1e6 < 1e-1:
                if verbose:
                    print("Problem calculating sensitivity for alert index {}".format(ind))
        except (IOError, ValueError, IndexError) as err:
            print(err)
    if with_disc:
        return sensitivities, discoveries
    else:
        return sensitivities

def ns_fits_contours(index, delta_t, smear=True, levs = [5., 25., 50., 75., 95.]):
    smear_str = 'smeared/' if smear else 'norm_prob/'
    levs = np.array(levs)
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/fits/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    with open(fs[0], 'r') as f:
        signal_trials = pickle.load(f)
    true_inj = np.array(signal_trials['true_ns'])
    ns_fit = np.array(signal_trials['ns_prior'])
    ninjs = np.unique(true_inj)
    if max(ninjs) < 10:
        print('Index {} has max {}'.format(index, max(ninjs)))
    contours = np.stack([np.percentile(ns_fit[true_inj==ninj], levs) for ninj in ninjs])
    return ninjs, contours.T

def ns_fits_contours_plot(index, delta_t, smear=True, levs=[5., 25., 50., 75., 95.],
                      show=False, col='navy green', custom_label = 'Median', ax=None,
                      xlabel=True, ylabel=True, legend=True):
    if ax is None:
        fig, ax = plt.subplots()
    ninj, fits = ns_fits_contours(index, delta_t, smear=smear, levs=levs)
    ax.plot(ninj, fits[2], label = custom_label, color = sns.xkcd_rgb[col])
    ax.fill_between(ninj, fits[0], fits[-1], alpha=0.3,
                    label='Central 90\%', color = sns.xkcd_rgb[col], lw=0)
    ax.fill_between(ninj, fits[1], fits[-2], alpha=0.5,
                    label='Central 50\%', color = sns.xkcd_rgb[col], lw=0)
    expectation = ninj
    exp_col = 'dark grey'
    ax.plot(ninj, expectation, ls = '--', color = sns.xkcd_rgb[exp_col])
    if legend:
        ax.legend(loc=4)
    if xlabel:
        ax.set_xlabel(r'$n_{\mathrm{inj}}$')
    ax.set_xlim(0., max(ninj))
    ax.set_ylim(0., 80)
    if ylabel:
        ax.set_ylabel(r'$\hat{n}_{s}$')
    if show:
        plt.show()

def fitting_bias_summary(delta_t, sigs=[2., 5., 10.], smear=True, containment=50.):
    bias = {sig: [] for sig in sigs}; spread = {sig: [] for sig in sigs};
    levs = [50.-containment / 2., 50., 50.+containment / 2.]
    for ind in range(249):
        try:
            ninjs, contours = ns_fits_contours(ind, delta_t, smear=smear, levs=levs)
        except:
            for sig in sigs:
                bias[sig].append(0.0)
                spread[sig].append(0.0)
            continue
        for sig in sigs:
            try:
                n_ind = np.argwhere(ninjs == sig)[0][0]
                bias[sig].append(contours[1][n_ind])
                spread[sig].append(contours[-1][n_ind] - contours[0][n_ind])
            except:
                bias[sig].append(0.0)
                spread[sig].append(0.0)
    return bias, spread

def background(index, delta_t, smear=True):
    smear_str = 'smeared/' if smear else 'norm_prob/'
    fs = glob('/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/bg/{}index_{}_*_time_{:.1f}.pkl'.format(smear_str, index, delta_t))
    with open(fs[0], 'r') as f:
        bg_trials = pickle.load(f)
    return bg_trials

def plot_zoom_from_map(skymap, ind, reso=1., cmap=None, draw_contour=True, ax=None, 
                      col_label= r'$\log_{10}$(prob.)'):
    s, header = hp.read_map(skymap_files[ind], h=True, verbose=False)
    header = {name: val for name, val in header}
    nside = hp.get_nside(s)
    area = np.count_nonzero(s < 64.2) * hp.nside2pixarea(nside) * 180.**2. / (np.pi**2.)
    reso *= int(np.sqrt(area))
    reso = np.max([reso, 1.])
    original_LLH = s
    ra = np.radians(header['RA'])
    dec = np.radians(header['DEC'])
    title = skymap_files[ind][l_ind:r_ind].replace('Run', 'Run ').replace('_', ', Event ')
    if cmap is None:
        pdf_palette = sns.color_palette("Blues", 500)
        cmap = mpl.colors.ListedColormap(pdf_palette)
    if np.count_nonzero(skymap > 0.0) > 1:
        max_color = np.max(skymap)
        min_color = 0.
    else:
        max_color =  -1.8 #5 #max(skymap)
        min_color = -5.  #0.
    #min_color = np.min([0., 2.*max_color])
    hp.gnomview(skymap, rot=(np.degrees(ra), np.degrees(dec), 0),
                    cmap=cmap,
                    max=max_color,
                    min=min_color,
                    reso=reso,
                    title=title,
                    notext=True,
                    cbar=False
                    #unit=r""
                    )

    plt.plot(4.95/3.*reso*np.radians([-1, 1, 1, -1, -1]), 4.95/3.*reso*np.radians([1, 1, -1, -1, 1]), color="k", ls="-", lw=3)
    hp.graticule(verbose=False)
    steady_sensitivity_fits.plot_labels(dec, ra, reso)
    con_nside = 256 if area < 5. else 128
    if draw_contour:
        contours = steady_sensitivity_fits.plot_contours(None, original_LLH, levels=[22.2, 64.2], nside = con_nside)
        for contour in np.array(contours).T:
            hp.projplot(contour[0],contour[1],linewidth=1.5,c='k')
    steady_sensitivity_fits.plot_color_bar(cmap = cmap, labels = [min_color, max_color], col_label = col_label)


def background_hotspot_map(ind, delta_t, smear=True):
    bg = background(ind, delta_t, smear=smear)
    msk = np.array(bg['ts_prior']) != 0.
    ra, dec = np.array(bg['ra'])[msk], np.array(bg['dec'])[msk]
    theta = np.pi/2. - dec
    inds = hp.ang2pix(256, theta, ra)
    ind_counts = np.unique(inds, return_counts=True)
    reco_hist = np.zeros(hp.nside2npix(256))
    reco_hist[ind_counts[0]] = ind_counts[1]
    plot_zoom_from_map(reco_hist, ind, draw_contour=False, col_label='Counts')

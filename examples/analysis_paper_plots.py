import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
import pickle
from matplotlib.lines import Line2D
from matplotlib.legend import Legend
from matplotlib.colors import ListedColormap

import scipy as sp
#from lmfit                import Model
from scipy.optimize       import curve_fit
from scipy.stats          import chi2
from scipy import ndimage

mpl.style.use('/home/apizzuto/Nova/scripts/novae_plots.mplstyle')
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['figure.dpi'] = 120
mpl.rcParams['font.family'] = 'Times New Roman'

#sns.palplot(sns.color_palette('colorblind'))
palette = ['#8da0cb', '#fc8d62', '#66c2a5']

def erfunc(x, a, b):
    return 0.5 + 0.5*sp.special.erf(a*x + b)

def chi2cdf(x,df1,loc,scale):
    func = chi2.cdf(x,df1,loc,scale)
    return func

def incomplete_gamma(x, a, scale):
    return sp.special.gammaincc( scale*x, a)

def poissoncdf(x, mu, loc):
    func = sp.stats.poisson.cdf(x, mu, loc)
    return func

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def n_to_flux(N, deltaT, sinDec, gamma, month=6):
    signal_trials = np.load('/data/user/apizzuto/fast_response_skylab/fast-response/fast_response/analysis_checks/sensitivity/nsignal_sinDec_{}_deltaT_{}_month_{}.npy'.format(sinDec, deltaT, month))
#     with open('/data/user/apizzuto/fast_response_skylab/alert_event_followup/sensitivity_ts_distributions/index_{}_time_{:.1f}.pkl'.format(index, delta_t), 'r') as f:
#         signal_trials = pickle.load(f)
    signal_trials = signal_trials[signal_trials['gamma'] == gamma]
    fl_per_one = np.mean(signal_trials['flux'] / signal_trials['mean_ninj'])
    return fl_per_one * N

def pass_vs_inj(deltaT, sinDec, gamma, threshold = 0.5, in_ns = True, month=6, with_err = True, trim=-1):
    with open('/data/user/apizzuto/fast_response_skylab/fast-response/fast_response/analysis_checks/bg/ps_sinDec_{}_deltaT_{}_month_{}.pkl'.format(sinDec, deltaT, month), 'r') as f:
        bg_trials = pickle.load(f)
    bg_trials = np.array([0.0]*bg_trials['n_zero'] + list(bg_trials['tsd']))
    signal_trials = np.load('/data/user/apizzuto/fast_response_skylab/fast-response/fast_response/analysis_checks/sensitivity/nsignal_sinDec_{}_deltaT_{}_month_{}.npy'.format(sinDec, deltaT, month))
    signal_trials = signal_trials[signal_trials['gamma'] == gamma]
    bg_thresh = np.percentile(bg_trials, threshold * 100.)
    signal_fluxes, signal_indices = np.unique(signal_trials['mean_ninj'], return_index=True)
    signal_indices = np.append(signal_indices, len(signal_trials))
    if trim != -1 and trim < 0:
        signal_indices = signal_indices[:trim]
        signal_fluxes = signal_fluxes[:trim]
    elif trim > 0:
        signal_indices = signal_indices[:trim + 1]
        signal_fluxes = signal_fluxes[:trim]
    passing = np.array([np.count_nonzero(signal_trials['TS'][li:ri] > bg_thresh) / float(ri - li) for li, ri in zip(signal_indices[:-1], signal_indices[1:])])
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
    
def sensitivity_curve(deltaT, sinDec, gamma, threshold = 0.5, in_ns = True, month=6, with_err = True, trim=-1, ax = None, p0 = None, fontsize = 16):
    signal_fluxes, passing, errs = pass_vs_inj(deltaT, sinDec, gamma, threshold=threshold, in_ns=in_ns, month=month, with_err=with_err, trim=trim)
    fits, plist = [], []
    try:
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, chi2cdf, p0=p0))
        plist.append(fits[-1]['pval'])
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, erfunc, p0=p0))
        plist.append(fits[-1]['pval'])
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, incomplete_gamma, p0=p0))
        plist.append(fits[-1]['pval'])
    except:
        pass
        #print("at least one fit failed")
    #Find best fit of the three, make it look different in plot
    plist = np.array(plist)
    best_fit_ind= np.argmax(plist)
    fits[best_fit_ind]['ls'] = '-'
    
    if ax==None:
        fig, ax = plt.subplots()
    
    for fit_dict in fits:
        ax.plot(fit_dict['xfit'], fit_dict['yfit'], 
                 label = r'{}: $\chi^2$ = {:.2f}, d.o.f. = {}'.format(fit_dict['name'], fit_dict['chi2'], fit_dict['dof']),
                ls = fit_dict['ls'])
        if fit_dict['ls'] == '-':
            ax.axhline(0.9, color = palette[-1], linewidth = 0.3, linestyle = '-.')
            ax.axvline(fit_dict['sens'], color = palette[-1], linewidth = 0.3, linestyle = '-.')
            ax.text(5, 0.8, r'Sens. = {:.2f}'.format(fit_dict['sens']))
    ax.errorbar(signal_fluxes, passing, yerr=errs, capsize = 3, linestyle='', marker = 's', markersize = 2)
    ax.legend(loc=4, fontsize = fontsize)
    
def calc_sensitivity(deltaT, sinDec, gamma, threshold = 0.5, in_ns = True, month=6, with_err = True, trim=-1, 
                     p0=None, conf_lev = 0.9):
    signal_fluxes, passing, errs = pass_vs_inj(deltaT, sinDec, gamma, threshold=threshold, in_ns=in_ns, month=month, with_err=with_err, trim=trim)
    fits, plist = [], []
    try:
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, chi2cdf, p0=p0, conf_lev=conf_lev))
        plist.append(fits[-1]['pval'])
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, erfunc, p0=p0, conf_lev=conf_lev))
        plist.append(fits[-1]['pval'])
        fits.append(sensitivity_fit(signal_fluxes, passing, errs, incomplete_gamma, p0=p0, conf_lev=conf_lev))
        plist.append(fits[-1]['pval'])
    except:
        pass
    #Find best fit of the three, make it look different in plot
    plist = np.array(plist)
    best_fit_ind= np.argmax(plist)
    return fits[best_fit_ind]
    
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

def pvals_for_signal(deltaT, sinDec, gamma, ns = 1, month=6, sigma = False):
    with open('/data/user/apizzuto/fast_response_skylab/fast-response/fast_response/analysis_checks/bg/ps_sinDec_{}_deltaT_{}_month_{}.pkl'.format(sinDec, deltaT, month), 'r') as f:
        bg_trials = pickle.load(f)
    bg_trials = np.array([0.0]*bg_trials['n_zero'] + list(bg_trials['tsd']))
    signal_trials = np.load('/data/user/apizzuto/fast_response_skylab/fast-response/fast_response/analysis_checks/sensitivity/nsignal_sinDec_{}_deltaT_{}_month_{}.npy'.format(sinDec, deltaT, month))
    signal_trials = signal_trials[signal_trials['gamma'] == gamma]
    signal_trials = signal_trials[signal_trials['n_inj'] == ns]
    #print(len(bg_trials['TS']))
    pvals = [100. - sp.stats.percentileofscore(bg_trials, ts, kind='strict') for ts in signal_trials['TS']]
    pvals = np.array(pvals)*0.01
    min_allowed_p = 1. / (float(len(bg_trials)) * 2.)
    pvals = np.where(pvals==0, min_allowed_p, pvals)
    if not sigma:
        return pvals
    else:
        return sp.stats.norm.ppf(1. - (pvals / 2.))

deltaTs = np.logspace(3., 7., 3)
sinDecs = np.linspace(-0.9, 0.9, 19)
gamma = 2.0

sensitivity = {}
discovery = {}
for jj, deltaT in enumerate(deltaTs):
    sensitivity[deltaT] = []
    discovery[deltaT] = []
    for ii, sinDec in enumerate(sinDecs):
        fit_dict = calc_sensitivity(deltaT, sinDec, gamma)
        flux_conv = n_to_flux(1, deltaT, sinDec, gamma, month=6)
        sensitivity[deltaT].append(fit_dict['sens']*flux_conv*deltaT*1e6)
        fit_dict = calc_sensitivity(deltaT, sinDec, gamma, threshold=1.0-0.0013, conf_lev=0.9)
        discovery[deltaT].append(fit_dict['sens']*flux_conv*deltaT*1e6)

sens_cols = palette
#disc_cols = [sns.xkcd_rgb['light navy blue'], sns.xkcd_rgb['lighter green'], sns.xkcd_rgb['light magenta']]
fig, axs = plt.subplots(ncols=3, nrows=1, sharey=True, figsize = (13, 4), dpi=200)
plt.subplots_adjust(wspace=0.08)
fig.set_facecolor('white')

for ii, delta_t in enumerate(deltaTs[:]):
    axs[ii].plot(sinDecs, sensitivity[delta_t], lw=3, ls='-', 
             color=sens_cols[ii], label = 'Sensitivity')
    axs[ii].plot(sinDecs,discovery[delta_t], lw=3, ls='--', 
             color=sens_cols[ii], label = r'$3\sigma$ discovery') # +'\n\t(50\% CL)')

    axs[ii].grid(which='both', alpha=0.2, zorder=1)
    axs[ii].set_yscale('log')
    axs[ii].legend(loc=1, ncol=1, frameon=False) #columnspacing=0.6, frameon=False)
    #plt.loglog()
    axs[ii].set_xlabel(r'$\sin \delta$')
    exp = np.log10(delta_t)
    labs = [r'$\Delta T = 10^{3}\,$s', r'$\Delta T = 10^{5}\,$s', r'$\Delta T = 10^{7}\,$s']
    axs[ii].set_title(labs[ii]) # + '{:.1e} s'.format(delta_t))
    if ii == 0:
        axs[ii].set_ylabel(r'$E^2 \frac{dN_{\nu_{\mu}+\bar{\nu}_{\mu}}}{dEdA} \Big|_{\mathrm{1 TeV}}$ (GeV cm$^{-2}$)')
    axs[ii].set_ylim(8e-3, 5e0)
    #plt.text(0.03, 1e0, '10 yr. time-integrated', color = sns.xkcd_rgb['light grey'])
    axs[ii].set_xlim(-1.05, 1.05)
    axs[ii].set_xticks(np.linspace(-1., 1., 5))
plt.savefig('figures/sensitivity_time_window_panel_plot.png', 
           dpi=200, bbox_inches ='tight')

deltaTs = np.logspace(1., 7.5, 14)
sinDecs = np.array([-0.5, 0.0, 0.5])
gamma = 2.0

sensitivity = {}
discovery = {}
for ii, sinDec in enumerate(sinDecs):
    sensitivity[sinDec] = []
    discovery[sinDec] = []
    for jj, deltaT in enumerate(deltaTs):
        try:
            fit_dict = calc_sensitivity(deltaT, sinDec, gamma)
            flux_conv = n_to_flux(1, deltaT, sinDec, gamma, month=6)
            sensitivity[sinDec].append(fit_dict['sens']*flux_conv*deltaT*1e6)
            fit_dict = calc_sensitivity(deltaT, sinDec, gamma, threshold=1.0-0.0013, conf_lev=0.9)
            discovery[sinDec].append(fit_dict['sens']*flux_conv*deltaT*1e6)
        except:
            print(deltaT)
            pass

sens_cols = palette
fig, ax = plt.subplots(dpi=200, figsize = (8, 4))
fig.set_facecolor('white')
custom_lines = [Line2D([0], [0], color=sns.xkcd_rgb['dark grey'], lw=3, ls = '-'),
                Line2D([0], [0], color=sns.xkcd_rgb['dark grey'], lw=3, ls = '--')]

for jj, sinDec in enumerate(sinDecs):
    try:
        msk = deltaTs != 316.22776601683796
        plt.plot(deltaTs[msk], sensitivity[sinDec], lw=3, ls='-', 
                 color=sens_cols[jj], label = r'$\sin \delta = $' + ' {:.1f}'.format(sinDec))
        plt.plot(deltaTs[msk], discovery[sinDec], lw=3, ls='--', 
                 color=sens_cols[jj]) # , label = r'$3\sigma$ discovery' +'\n\t(90\% CL)')
    except:
        plt.plot(deltaTs, sensitivity[sinDec], lw=3, ls='-', 
                 color=sens_cols[jj], label = r'$\sin \delta = $' + ' {:.1f}'.format(sinDec))
        plt.plot(deltaTs, discovery[sinDec], lw=3, ls='--', 
                 color=sens_cols[jj]) # , label = r'$3\sigma$ discovery' +'\n\t(90\% CL)')

plt.grid(which='both', alpha=0.2, zorder=1)
plt.yscale('log')
plt.legend(loc=2, ncol=3, frameon=False, ) #columnspacing=0.6, frameon=False)
plt.xlabel(r'$\Delta T$ (s)')
plt.ylabel(r'$E^2 \frac{dN_{\nu_{\mu}+\bar{\nu}_{\mu}}}{dEdA} \Big|_{\mathrm{1 TeV}}$ (GeV cm$^{-2}$)')
plt.ylim(8e-3, 5e0)
plt.xscale('log')

leg = Legend(ax, custom_lines, ['Sensitivity', r'$3\sigma$ discovery'], # +' (90\% CL)'],
             loc='lower center', frameon=False, ncol = 2)
ax.add_artist(leg)
plt.savefig('figures/sensitivity_declination_panel_plot_90CL.png', 
           dpi=200, bbox_inches ='tight')

fig, ax = plt.subplots(dpi=200, figsize = (8, 4))
fig.set_facecolor('white')
custom_lines = [Line2D([0], [0], color=sns.xkcd_rgb['dark grey'], lw=3, ls = '-'),
                Line2D([0], [0], color=sns.xkcd_rgb['dark grey'], lw=3, ls = '--')]

for jj, sinDec in enumerate(sinDecs):
    try:
        msk = deltaTs != 316.22776601683796
        plt.plot(deltaTs[msk], sensitivity[sinDec] / deltaTs[msk], lw=3, ls='-', 
                 color=sens_cols[jj], label = r'$\sin \delta = $' + ' {:.1f}'.format(sinDec))
        plt.plot(deltaTs[msk], discovery[sinDec] / deltaTs[msk], lw=3, ls='--', 
                 color=sens_cols[jj]) # , label = r'$3\sigma$ discovery' +'\n\t(90\% CL)')
    except:
        plt.plot(deltaTs, sensitivity[sinDec] / deltaTs, lw=3, ls='-', 
                 color=sens_cols[jj], label = r'$\sin \delta = $' + ' {:.1f}'.format(sinDec))
        plt.plot(deltaTs, discovery[sinDec] / deltaTs, lw=3, ls='--', 
                 color=sens_cols[jj]) # , label = r'$3\sigma$ discovery' +'\n\t(90\% CL)')

plt.grid(which='both', alpha=0.2, zorder=1)
plt.yscale('log')
plt.legend(loc=1, ncol=1, frameon=False, ) #columnspacing=0.6, frameon=False)
plt.xlabel(r'$\Delta T$ (s)')
plt.ylabel(r'$E^2 \frac{dN_{\nu_{\mu}+\bar{\nu}_{\mu}}}{dEdA} \Big|_{\mathrm{1 TeV}}$ (GeV cm$^{-2}$ s$^{-1}$)')
#plt.ylim(8e-3, 5e0)
plt.xscale('log')

leg = Legend(ax, custom_lines, ['Sensitivity', r'$3\sigma$ discovery'], # +' (90\% CL)'],
             loc='lower center', frameon=False, ncol = 2)
ax.add_artist(leg)
plt.xlim(1e2, 5e5)
plt.ylim(3e-8, 1e-2)
fig.set_facecolor('w')
plt.savefig('figures/sensitivity_vs_time_in_flux.png', 
           dpi=200, bbox_inches ='tight')

deltaTs = np.logspace(1., 7.5, 14)
sinDecs = np.array([-0.5, 0.0, 0.5])
gamma = 2.0

sensitivity = {}
discovery = {}
for ii, sinDec in enumerate(sinDecs):
    sensitivity[sinDec] = []
    discovery[sinDec] = []
    for jj, deltaT in enumerate(deltaTs):
        try:
            fit_dict = calc_sensitivity(deltaT, sinDec, gamma)
            flux_conv = n_to_flux(1, deltaT, sinDec, gamma, month=6)
            sensitivity[sinDec].append(fit_dict['sens']*flux_conv*deltaT*1e6)
            fit_dict = calc_sensitivity(deltaT, sinDec, gamma, threshold=1.0-0.0013, conf_lev=0.5)
            discovery[sinDec].append(fit_dict['sens']*flux_conv*deltaT*1e6)
        except:
            print(deltaT)
            pass

msk = deltaTs != 316227.7660168379
fig, ax = plt.subplots(dpi=200, figsize = (8, 4))
fig.set_facecolor('white')
custom_lines = [Line2D([0], [0], color=sns.xkcd_rgb['dark grey'], lw=3, ls = '-'),
                Line2D([0], [0], color=sns.xkcd_rgb['dark grey'], lw=3, ls = '--')]

for jj, sinDec in enumerate(sinDecs):
    try:
        msk = deltaTs != 316.22776601683796
        plt.plot(deltaTs[msk], sensitivity[sinDec], lw=3, ls='-', 
                 color=sens_cols[jj], label = r'$\sin \delta = $' + ' {:.1f}'.format(sinDec))
        plt.plot(deltaTs[msk], discovery[sinDec], lw=3, ls='--', 
                 color=sens_cols[jj]) # , label = r'$3\sigma$ discovery' +'\n\t(90\% CL)')
    except:
        plt.plot(deltaTs, sensitivity[sinDec], lw=3, ls='-', 
                 color=sens_cols[jj], label = r'$\sin \delta = $' + ' {:.1f}'.format(sinDec))
        plt.plot(deltaTs, discovery[sinDec], lw=3, ls='--', 
                 color=sens_cols[jj]) # , label = r'$3\sigma$ discovery' +'\n\t(90\% CL)')

plt.grid(which='both', alpha=0.2, zorder=1)
plt.yscale('log')
plt.legend(loc=2, ncol=3, frameon=False, ) #columnspacing=0.6, frameon=False)
plt.xlabel(r'$\Delta T$ (s)')
plt.ylabel(r'$E^2 \frac{dN_{\nu_{\mu}+\bar{\nu}_{\mu}}}{dEdA} \Big|_{\mathrm{1 TeV}}$ (GeV cm$^{-2}$)')
plt.ylim(2e-3, 5e0)
plt.xscale('log')

leg = Legend(ax, custom_lines, ['Sensitivity', r'$3\sigma$ discovery'], # +' (50\% CL)'],
             loc='lower center', frameon=False, ncol = 2)
ax.add_artist(leg)
plt.savefig('figures/sensitivity_declination_panel_plot_50CL.png', 
           dpi=200, bbox_inches ='tight')

fig = plt.figure(dpi=200)
plt.subplots_adjust(wspace=0.05)
fig.set_facecolor('w')
bin_width = np.sqrt(10.)

sens_cols = palette

delta_t = 1e3
for ii, (ind, dec, mark) in enumerate(zip([1, 2, 5], [-30., 0., 30.], ['^', '8', 's'])):
    ndays = delta_t / 86400.
    gfu = np.load('/data/user/apizzuto/fast_response_skylab/dump/differential_sensitivity_deltaT_{:.2e}_sinDec_{:.2f}.pkl'.format(ndays, np.sin(np.radians(dec))),
               allow_pickle=True)
    ens = np.array(gfu['low_energies'])*np.sqrt(10.)
    sens = np.array(gfu['sensitivity']) * 1e3 #1e6 for E_0^2, 1e-3 for GeV to TeV
    sens *= 1e3 * delta_t #Back to GeV integrate over time
    plt.errorbar(ens, sens,  
        xerr=[ens-ens/bin_width, ens*bin_width-ens],
        marker=mark, ls=':', label=r'$\sin \delta = $' +'{:.1f}'.format(np.sin(np.radians(dec))),
                color=sens_cols[ii])
plt.xscale('log')
plt.yscale('log')
plt.ylim(3e-2, 1e4)
plt.xlim(3e1, 2e9)
plt.legend(loc='upper left', fontsize=12, ncol=1)
plt.xlabel(r'$E_{\nu}$ (GeV)')
plt.text(5e1, 6e1, r'$\Delta T = 10^3$ s') #+'{:.1e} s'.format(delta_t))
#if ii == 0:
plt.ylabel(r'$E^2 \frac{dN_{\nu_{\mu}+\bar{\nu}_{\mu}}}{dEdA} \Big|_{\mathrm{1 TeV}}$ (GeV cm$^{-2}$)')
        
#fig.suptitle(r'Sensitivity to $\Delta t =$' + '{:.1e} s transients'.format(delta_t), y=0.95)
plt.savefig('figures/differential_sensitivity_delta_t_{:.1e}.png'.format(delta_t),
           bbox_inches='tight')

deltaTs = np.logspace(1.5, 6.5, 11)
sinDecs = np.linspace(-0.9, 0.9, 37)

for gamma in [2.0, 2.5, 3.0]:

    med_sigs = []
    for ii, sinDec in enumerate(sinDecs):
        print(ii),
        med_sigs.append([])
        for jj, deltaT in enumerate(deltaTs):
            #if jj%4==0:
                #print jj
            try:
                significances = pvals_for_signal(deltaT, sinDec, gamma, ns = 1, sigma = True, month = 6)
                med_sigs[-1].append(np.median(significances))
                srbs = significances
            except IOError, e:
                med_sigs[-1].append(np.median(srbs))
                print(e)

    sigma = 1.5 # this depends on how noisy your data is, play with it!
    data = ndimage.filters.gaussian_filter(med_sigs, sigma)

    fig, ax = plt.subplots()

    cmap = ListedColormap(sns.cubehelix_palette(12, start=.5, rot=-.75, reverse=True))
    X, Y = np.meshgrid(np.linspace(1.5, 6.5, 11), sinDecs)

    cs = ax.contour(X,Y, data, np.linspace(0., 5.5, 12), cmap=cmap)
    csf = ax.contourf(X,Y, data, np.linspace(0., 5.5, 12), cmap=cmap)
    ax.clabel(cs, cs.levels[::2], inline=True, fontsize=16, colors = 'k', 
            inline_spacing=12, fmt = r'%1.1f $\sigma$', rightside_up=True)
    cbar = plt.colorbar(csf, ticks=np.linspace(1., 5., 5))
    cbar.set_label(r'Median significance ($\sigma$)', fontsize = 18)
    cbar.ax.tick_params(axis='y', direction='out')
    plt.xlabel(r'$\log_{10} \Big(\Delta T / \rm{s} \Big)$', fontsize = 18)
    plt.title('One Injected Event, $\gamma = ${}'.format(gamma), fontsize = 20)
    plt.ylabel(r'$\sin(\delta)$', fontsize = 18)
    plt.xlim(1.9, 6.4)
    plt.savefig('figures/FRA_significance_gamma_{}_one_event.png'.format(gamma), dpi=200, bbox_inches='tight')


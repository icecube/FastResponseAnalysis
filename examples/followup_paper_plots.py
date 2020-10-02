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

mpl.style.use('/home/apizzuto/Nova/scripts/novae_plots.mplstyle')

#### BEGIN WITH 190114C LIMITS ####

gev_to_erg = 0.00160218
icecube_delta_t = 3750.
icecube_lim = 3.5e-01 * gev_to_erg / icecube_delta_t
icecube_es = np.array([8e4, 2e7]) * 1e9
mid_icecube_e = 10.**np.mean(np.log10(icecube_es))

antares_delta_t = 1600.
antares_es = np.array([1e4, 2e7]) * 1e9
antares_lim = 1.6 * gev_to_erg / antares_delta_t
mid_antares_e = 10.**np.mean(np.log10(antares_es))

def load_band(e_range, time_window):
    low_edge = np.genfromtxt('/data/user/apizzuto/fast_response_skylab/dump/grb190114c_data/time_window_{}_{}_bottom.csv'.format(time_window, e_range),
             delimiter=', ')
    high_edge = np.genfromtxt('/data/user/apizzuto/fast_response_skylab/dump/grb190114c_data/time_window_{}_{}_top.csv'.format(time_window, e_range),
             delimiter=', ')
    low_spl = sp.interpolate.UnivariateSpline(np.log10(low_edge.T[0]), np.log10(low_edge.T[1]))
    high_spl = sp.interpolate.UnivariateSpline(np.log10(high_edge.T[0]), np.log10(high_edge.T[1]))
    low_e = np.max([low_edge.T[0].min(), high_edge.T[0].min()])
    high_e = np.min([low_edge.T[0].max(), high_edge.T[0].max()])
    es = np.linspace(low_e, high_e, 500)
    return es, 10.**low_spl(np.log10(es)), 10.**high_spl(np.log10(es))

def load_bands(time_window):
    low = load_band('low', time_window)
    med = load_band('med', time_window)
    high = load_band('hi', time_window)
    return low, med, high

def load_points(time_window):
    data = np.genfromtxt('/data/user/apizzuto/fast_response_skylab/dump/grb190114c_data/time_window_{}_data.csv'.format(time_window),
             delimiter=', ')
    xerr0 = np.genfromtxt('/data/user/apizzuto/fast_response_skylab/dump/grb190114c_data/time_window_{}_left.csv'.format(time_window),
             delimiter=', ').T[0]
    xerr1 = np.genfromtxt('/data/user/apizzuto/fast_response_skylab/dump/grb190114c_data/time_window_{}_right.csv'.format(time_window),
             delimiter=', ').T[0] 
    yerr0 = np.genfromtxt('/data/user/apizzuto/fast_response_skylab/dump/grb190114c_data/time_window_{}_down.csv'.format(time_window),
             delimiter=', ').T[1]
    yerr1 = np.genfromtxt('/data/user/apizzuto/fast_response_skylab/dump/grb190114c_data/time_window_{}_up.csv'.format(time_window),
             delimiter=', ').T[1]
    return data, [data.T[0] - xerr0, xerr1 - data.T[0]], [data.T[1] - yerr0, yerr1 - data.T[1]]

fig = plt.figure(dpi=200)
bin_width = np.sqrt(10.)
fig.set_facecolor('w')

#plt.errorbar(ens, sens,  xerr=[ens-ens/bin_width, ens*bin_width-ens],
#                marker='v', ls=':', color=sns.xkcd_rgb['dark sky blue'])
plt.errorbar([mid_icecube_e], [icecube_lim], xerr=[np.array([mid_icecube_e]) - icecube_es[0],
                    -1.*np.array([mid_icecube_e]) + icecube_es[1]],
                    marker='v', color=sns.xkcd_rgb['dark magenta'], 
                     label = 'IceCube U.L. (-150 s, 3600 s)')

plt.errorbar([mid_antares_e], [antares_lim], xerr=[np.array([mid_antares_e]) - antares_es[0],
                    -1.*np.array([mid_antares_e]) + antares_es[1]],
                    marker='v', color='grey', label = 'ANTARES U.L. (-350 s, 1250 s)')

cols = [sns.xkcd_rgb['dark navy blue'], sns.xkcd_rgb['dark sky blue']] #dark sky blue, deep sky blue
time_strs = ['Swift, $\it{Fermi}$, MAGIC (68 s, 110 s)', 'Swift, $\it{Fermi}$, MAGIC (110 s, 180 s)']
for ii, time_window in enumerate([1, 2]):
    ds, xerr, yerr = load_points(time_window)
    plt.errorbar(ds.T[0], ds.T[1], xerr=xerr, yerr=yerr, ls='', 
                 marker='o', markersize=3, color = cols[ii], label=time_strs[ii])
    bands = load_bands(time_window)
    for band in bands:
        plt.fill_between(band[0], band[1], band[2], alpha=0.5, 
                         color = cols[ii], lw=0.)
#ds, xerr, yerr = load_points(2)
#plt.errorbar(ds.T[0], ds.T[1], xerr=xerr, yerr=yerr, ls='', marker='o', markersize=3)

plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-9, 5e-6)
plt.xlim(1e3, 1e17)
plt.legend(loc='upper left', fontsize=12, ncol=1, frameon=False)
plt.xlabel(r'$E$ (eV)')
plt.xticks(np.logspace(3, 15, 5))
#if ii == 0:
plt.ylabel(r'Flux (erg cm$^{-2}$ s$^{-1}$)')
plt.savefig('figures/GRB_190114C_upper_limit.png', bbox_inches='tight', dpi=200)


#### NOW PKS Modeling ####
gev_to_erg = 0.00160218
icecube_delta_t = 4.*86400. + (20.*3600)
icecube_lim =5.9e-01 * gev_to_erg / icecube_delta_t
icecube_es = np.array([8e4, 2e7]) * 1e9
mid_icecube_e = 10.**np.mean(np.log10(icecube_es))

icecube_sens =3.5e-2 * gev_to_erg / icecube_delta_t
icecube_sens_es = np.array([1e3, 5e6]) * 1e9
mid_icecube_sens_e = 10.**np.mean(np.log10(icecube_sens_es))

hz_to_erg = 4.136e-15

def load_archival():
    data = np.genfromtxt('/data/user/apizzuto/fast_response_skylab/dump/archival_pks0346.csv', delimiter=', ')
    inds = {'ATCA': (0, 3, 'o'), 'Planck': (4, 9, 's'), 'WISE': (10, 13, '^'), 'UVOT': (14, 26, '8'), 
            'XRT': (27, 29, 'p'), 'LAT': (30, 33, 'd')}
    up_lims = (34, 37)
    ens = data.T[0] * hz_to_erg
    fls = data.T[1]
    return ens, fls, inds, up_lims

def load_flare():
    data = np.genfromtxt('/data/user/apizzuto/fast_response_skylab/dump/flare_pks0346.csv', delimiter=', ')
    inds = {'UVOT': (0, 5, '8'), 'XRT': (6, 12, 'p'), 'LAT': (13, 15, 'd')}
    up_lims = (16, 18)
    ens = data.T[0] * hz_to_erg
    fls = data.T[1]
    return ens, fls, inds, up_lims

def load_cta():
    data = np.genfromtxt('/data/user/apizzuto/fast_response_skylab/dump/cta_sens.csv', delimiter=', ')
    ens = data.T[0] * hz_to_erg
    fls = data.T[1]
    return ens, fls

fig = plt.figure(dpi=200)
fig.set_facecolor('w')

ens, fls, inds, up_lims = load_archival()
for key in ['ATCA', 'Planck', 'WISE', 'UVOT', 'XRT', 'LAT']:
    plt.plot(ens[inds[key][0]:inds[key][1]+1], fls[inds[key][0]:inds[key][1]+1], ls='', marker=inds[key][2], 
             label=key, color=sns.xkcd_rgb['grey'], markeredgecolor='k')
plt.errorbar(ens[up_lims[0]:up_lims[1]+1], fls[up_lims[0]:up_lims[1]+1], ls='', marker=inds['LAT'][2], 
             yerr=np.array(fls[up_lims[0]:up_lims[1]+1])*0.3, uplims=True,
             color=sns.xkcd_rgb['grey'], markeredgecolor='k')
    
ens, fls, inds, up_lims = load_flare()
for key in ['UVOT', 'XRT', 'LAT']:
    plt.plot(ens[inds[key][0]:inds[key][1]+1], fls[inds[key][0]:inds[key][1]+1], ls='', marker=inds[key][2],
            color = sns.xkcd_rgb['dark sky blue'], markeredgecolor='k', markeredgewidth=1.)
plt.errorbar(ens[up_lims[0]:up_lims[1]+1], fls[up_lims[0]:up_lims[1]+1], ls='', marker=inds['LAT'][2], 
             yerr=np.array(fls[up_lims[0]:up_lims[1]+1])*0.3, uplims=True,
             color=sns.xkcd_rgb['dark sky blue'], markeredgecolor='k') 

plt.errorbar([mid_icecube_e], [icecube_lim], xerr=[np.array([mid_icecube_e]) - icecube_es[0],
                    -1.*np.array([mid_icecube_e]) + icecube_es[1]],
                    marker='v', color=sns.xkcd_rgb['dark magenta'], 
                     label = 'IceCube U.L.')

plt.errorbar([mid_icecube_sens_e], [icecube_sens], xerr=[np.array([mid_icecube_sens_e]) - icecube_sens_es[0],
                    -1.*np.array([mid_icecube_sens_e]) + icecube_sens_es[1]],
                    marker='v', color=sns.xkcd_rgb['grey'], markersize=0.,
                     label = 'IceCube sens.,\n' + r'$\delta=0.0^{\circ}$')

ens, fls = load_cta()
plt.plot(ens, fls, color = 'k', label='CTA (5 h)', ls = '--')


plt.loglog()
plt.xlabel(r'$E$ (eV)')
plt.xticks(np.logspace(-3, 15, 7))
plt.ylabel(r'Flux (erg cm$^{-2}$ s$^{-1}$)')
plt.ylim(6e-15, 4e-8)
plt.legend(loc='upper left', fontsize=12, frameon=True, ncol=2, columnspacing=0.5)
plt.savefig('figures/PKS_0346_upper_limit.png', bbox_inches='tight', dpi=200)


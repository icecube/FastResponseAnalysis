'''Using an analysis of one galactic FRB, calculate
limits on a population of FRBs assuming standard
candle luminosities. Note: this script requires
a functioning version of flarestack'''

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from flarestack import EnergyPDF
from flarestack.cosmo import get_rate, get_diffuse_flux, calculate_transient_cosmology, get_diffuse_flux_contour
from astropy import units as u
import numpy as np
import logging
import seaborn as sns
import copy
logging.getLogger().setLevel("INFO")
#logging.getLogger().setLevel("ERROR")
import sys
#sys.path.append('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/')
from fast_response.FastResponseAnalysis import FastResponseAnalysis
#mpl.rcParams['font.family'] = 'sans.serif'
mpl.style.use('/home/apizzuto/Nova/python3/scripts/novae_plots_nb.mplstyle')

##############################################
# First, calculate limit for an E^-2.5 limit #
##############################################

source = {'Location': "293.74, 21.89", 
         'Start Time': '2020-04-27 18:00:00.0',
         'Stop Time': '2020-04-28 18:00:00.0', 
         'Name': 'SGR 1935+2154 limit calculation', 
         #'save': False
         }
f = FastResponseAnalysis(source['Location'], source['Start Time'], source['Stop Time'], **source)
f.unblind_TS()
#f.calc_pvalue(ntrials=500)
up_lims = {}
for gamma in [2.0, 2.5]:
    f.initialize_injector(gamma=gamma)
    up_lims[gamma] = f.upper_limit()

for spec, lim in up_lims.items():
    print(f"Limit at 1 TeV: E^2dN/dE = {lim*86400*1e6:.2e} GeV per cm^2 for an E^-{spec} spectrum")

central_energies = {2.0: (8e2, 9e5),
                   2.5: (224.26383143837145, 82450.66191768281)}

##################################################
# Use limits to integrate over source population #
##################################################

frb_rate = get_rate("frb")
print(f"Local FRB rate is {frb_rate(0.0):.2g}")

dist = 16 * u.kpc

norm_energy = 1.*u.GeV

lim_dicts = {2.0: {'gamma': 2.0, 'flux_norm_lim': up_lims[2.0] * 86400. * 1000.**2.0 * (u.GeV / u.cm**2) / (norm_energy)**2.,
                  'central_ens': central_energies[2.0]},
            2.5: {'gamma': 2.5, 'flux_norm_lim': up_lims[2.5] * 86400. * 1000.**2.5 / (u.GeV * u.cm**2.),
                   'central_ens': central_energies[2.5]}}

e_pdf_dicts = {}

for spec, lim_dict in lim_dicts.items():
    spectrum_gamma = spec
    # dN/dE
    #atel_flux_norm_lim = 5.2 * 10**-2. * (u.GeV / u.cm**2) / (norm_energy)**2.
    e_pdf_dict = {
        "energy_pdf_name": "power_law",
        "gamma": spectrum_gamma,
        "e_min_gev": lim_dict['central_ens'][0],
        "e_max_gev": lim_dict['central_ens'][1],
    }

    epdf = EnergyPDF.create(e_pdf_dict)
    e_pdf_dicts[spec] = e_pdf_dict
    
    e_lim = (lim_dict['flux_norm_lim'] * epdf.fluence_integral() * norm_energy**2 * 4 * np.pi * dist.to("cm")**2.).to("erg")
    print(f"Muon Neutrino Energy limit for SGR 1935+2154 is {e_lim:.2g} between {epdf.e_min:.2g} GeV and {epdf.e_max:.2g} GeV")

non_standard_candle_fac = 5e2
scaled_lim_dicts = copy.deepcopy(lim_dicts[2.5])
scaled_lim_dicts['flux_norm_lim'] *= non_standard_candle_fac

limits = [
    ("Standard Candle, $E^{-2}$", lim_dicts[2.0]['flux_norm_lim'], e_pdf_dicts[2.0]),
    ("SGR 1935+2154-like FRBs, (Standard Candle)", lim_dicts[2.5]['flux_norm_lim'], e_pdf_dicts[2.5]),
    ("SGR 1935+2154-like FRBs, " + r'($\mathcal{E}_{\nu}\propto \mathcal{E}_{\mathrm{radio}}$)', scaled_lim_dicts['flux_norm_lim'], e_pdf_dicts[2.5])
]

fit = "joint_15"

fig = plt.figure(dpi=200)
ax = plt.subplot(111)
fig.set_facecolor('w')

best_fit, upper_butterfly, lower_butterfly, e_range = get_diffuse_flux_contour(fit=fit)
plt.plot(e_range, best_fit(e_range) * e_range**2, label="IceCube Diffuse Flux",
            color = sns.xkcd_rgb['dark navy blue'], zorder=4)
plt.fill_between(e_range, upper_butterfly(e_range)* e_range**2, 
            lower_butterfly(e_range)* e_range**2, alpha=0.3, 
            color = sns.xkcd_rgb['dark navy blue'], linewidth=0.0, zorder=3)

def frb_rate_high(z):
    return frb_rate(z) * (8.78+7.23) / 7.23
def frb_rate_low(z):
    return frb_rate(z) * (7.23-6.13) / 7.23

cols = [sns.xkcd_rgb['dark sky blue'], sns.xkcd_rgb['steel grey']]
counter = 0

for label, mean_flux_norm_lim, e_pdf_dict in limits[1:]:
    epdf = EnergyPDF.create(e_pdf_dict)
    lim_e_pdf_dict = dict(e_pdf_dict)
    lim_e_pdf_dict["nu_flux_at_1_gev"] = mean_flux_norm_lim * 4 * np.pi * dist**2.

    integrated_nu_flux_1_gev = calculate_transient_cosmology(
        lim_e_pdf_dict, frb_rate, "frb_limit", zmax=8.0, diffuse_fit=fit,
    )
    
    ens = np.logspace(np.log10(epdf.e_min), np.log10(epdf.e_max), 3.)
    fluxes = integrated_nu_flux_1_gev.value*np.power(ens, -1.*e_pdf_dict['gamma'])*ens**2.
    
    low_rate_integrated = calculate_transient_cosmology(
        lim_e_pdf_dict, frb_rate_low, "frb_limit", zmax=8.0, diffuse_fit=fit,
    )
    low_rate_fl = low_rate_integrated.value*np.power(ens, -1.*e_pdf_dict['gamma'])*ens**2.
    high_rate_integrated = calculate_transient_cosmology(
        lim_e_pdf_dict, frb_rate_high, "frb_limit", zmax=8.0, diffuse_fit=fit,
    )
    high_rate_fl = high_rate_integrated.value*np.power(ens, -1.*e_pdf_dict['gamma'])*ens**2.
    
    plt.errorbar(ens, fluxes, yerr=0.25*fluxes, uplims=True, 
                 label=label, color=cols[counter], zorder=2)
    plt.fill_between(ens, low_rate_fl, high_rate_fl, alpha=0.5, linewidth=0.,
                    color=cols[counter], zorder=3)

plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"$E_{\nu}$ (GeV)")
plt.ylabel(r"$E_{\nu}^{2} \frac{dN}{dE}$ (GeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$)")
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.4),
#          ncol=2, fancybox=True, shadow=True, fontsize=12)
ax.legend(loc=1, ncol=1, frameon=True, fontsize=10.5)
plt.xlim(1e2, 5e6)
plt.ylim(7e-12, 3e-5)
plt.savefig('figures/SGR1935_2154_FRB_standard_candle_lims.png', bbox_inches='tight', dpi=200)

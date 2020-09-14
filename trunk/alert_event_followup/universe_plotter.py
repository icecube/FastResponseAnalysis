import numpy as np
from glob import glob
import pandas as pd
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
import healpy as hp
import scipy.stats as st
import scipy as sp
import pickle
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
mpl.style.use('/home/apizzuto/Nova/scripts/novae_plots.mplstyle')

#skymap_files = glob('/data/ana/realtime/alert_catalog_v2/2yr_prelim/fits_files/Run*.fits.gz')
skymap_files = glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run1*.fits.gz')
energy_density = {'transient': {'HB2006SFR': 4.8038e51, 
                        'MD2014SFR': 6.196e51,
                        'NoEvolution': 1.8364e52},
                'steady': {'HB2006SFR': 7.6735e43, 
                        'MD2014SFR': 9.897e+43,
                        'NoEvolution': 2.93335e44}}

class UniversePlotter():
    r'''
    Tool to make some helpful plots

    Parameters:
    -----------
        - delta_t (float): Analysis timescale
        - data_years (float): Years of alert events
        - lumi (str): Luminosity function (standard candle of lognormal)
        - evol (str): Evolution model (ie Hopkins and Beacom, Madau and Dickinson)
    '''
    def __init__(self, delta_t, data_years, lumi, evol, **kwargs):
        self.delta_t = delta_t
        self.transient = True if self.delta_t is not None else False
        self.data_years = data_years
        self.lumi = lumi
        self.evol = evol
        self.background_median_ts = None
        self.background_median_p = None
        self.smeared = kwargs.pop('smeared', True)
        self.ts_fills = [self.background_median_ts, 50.] if self.transient else [self.background_median_ts, 10.]
        self.lumi_str = {'SC': 'Standard Candle', 'LG': 'Log Normal'}[self.lumi]
        self.evol_str = {'HB2006SFR': 'Hopkins and Beacom 2006 SFR',
                            'MD2014SFR': 'Madau and Dickinson 2014 CSFH'}[self.evol]
        self.ts_path = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/ts_distributions/'
        self.steady_str = '_delta_t_{:.2e}'.format(self.delta_t) if self.transient else '_steady'
        self.time_window_per_year = (365. * 86400.) / (self.delta_t) if self.transient else 1.
        key = 'transient' if self.transient else 'steady'
        self.energy_density = energy_density[key][evol]
        self.no_evol_energy_density = energy_density[key]['NoEvolution']
        self.evol_lumi_str = 'evol_{}_lumi_{}'.format(self.evol, self.lumi)
        self.densities = np.logspace(-11., -6., 21)
        self.luminosities = np.logspace(49, 62, 53) if self.transient else np.logspace(49., 56., 29)
        if self.transient:
            low_energy_msk = self.luminosities * self.delta_t * self.densities[0] < self.energy_density * 100.
            high_energy_msk = self.luminosities * self.delta_t * self.densities[-1] > self.energy_density * 1e-4
            self.luminosities = self.luminosities[low_energy_msk * high_energy_msk]
        self.seconds_per_year = 365.*86400.
        self.med_TS = None
        self.med_p = None
        self.get_labels()
        self.sigs = None

    def two_dim_sensitivity_plot_ts(self, compare=False, log_ts=False, in_ts=True,
                        ts_vs_p=False):
        r'''Two dimensional contour plot to show the sensitivity of the analysis
        in the luminosity-density plane, highlighting the energy requirements
        for the diffuse flux

        Parameters
        ----------
            - compare (bool): include constraints from previous analyses
            - log_ts (bool): contour colors from log or linear ts values
            - in_ts (bool): TS value or binomial-p value
            - ts_vs_p (bool): compare TS and binomial-p value construction
        '''
        if in_ts or ts_vs_p:
            if self.background_median_ts is None:
                self.get_overall_background_ts()
            if self.med_TS is None:
                self.get_med_TS()
        if (not in_ts) or ts_vs_p:
            if self.background_median_p is None:
                self.get_overall_background_p()
            if self.med_p is None:
                self.get_med_p()
        fig, ax = plt.subplots(figsize=(8,5), dpi=200)
        fig.set_facecolor('w')
        X, Y = np.meshgrid(np.log10(self.densities), np.log10(self.plot_lumis))
        if in_ts or ts_vs_p:
            plot_vals = self.med_TS if not log_ts else np.log10(self.med_TS) 
        else:
            plot_vals = self.med_p if not log_ts else np.log10(self.med_p) 
        if in_ts:
            if log_ts:
                levels = np.linspace(np.log10(self.background_median_ts), np.log10(self.background_median_ts) + 2., 11)
            else:
                levels = np.linspace(self.background_median_ts, self.background_median_ts * 10., 11)
        else:
            if log_ts:
                levels = np.linspace(-8., np.log10(self.background_median_p), 11)
            else:
                levels = np.linspace(0.0, self.background_median_p, 11)
        cmap = self.cmap if in_ts else ListedColormap(self.cmap.colors[::-1])
        extend = 'max' if in_ts else 'min'
        cs = ax.contour(X, Y, plot_vals, cmap=cmap, levels=levels, 
                        #vmin=-0.5, 
                        extend=extend)
        csf = ax.contourf(X, Y, plot_vals, cmap=cmap, levels=levels, 
                        #vmin=-0.5, 
                        extend=extend)
        cbar = plt.colorbar(csf) 
        if in_ts or ts_vs_p:
            cbar_lab = r'Median Stacked TS' if not log_ts else r'$\log_{10}($Median Stacked TS$)$'
        else:
            cbar_lab = r'Median binom. p' if not log_ts else r'$\log_{10}($Median binom. p$)$'
        cbar.set_label(cbar_lab, fontsize = 18)
        cbar.ax.tick_params(axis='y', direction='out')
        if in_ts or ts_vs_p:
            cs_ts = ax.contour(X, Y, self.lower_10 - self.background_median_ts, colors=['k'], 
                            levels=[0.0], linewidths=2., zorder=10)
        if (not in_ts) or ts_vs_p:
            linestyle = 'solid' if ts_vs_p else 'dashed'
            cs_ts = ax.contour(X, Y, self.background_median_p - self.lower_10_p, colors=['k'], 
                            levels=[0.0], linewidths=2., linestyles=linestyle)
        xs = np.logspace(-11., -6., 1000)
        ys_max = self.no_evol_energy_density / xs / self.seconds_per_year if self.transient else self.no_evol_energy_density / xs
        ys_min = self.energy_density / xs / self.seconds_per_year if self.transient else self.energy_density / xs
        plt.fill_between(np.log10(xs), np.log10(ys_min), np.log10(ys_max), 
                color = 'm', alpha = 0.3, lw=0.0, zorder=10)
        if compare:
            comp_rho, comp_en, comp_str = self.compare_other_analyses()
            plt.plot(comp_rho, comp_en, color = 'gray', lw=2., zorder=5)
        plt.text(-9, 54.1, 'Diffuse', color = 'm', rotation=-28, fontsize=18)
        #plt.text(-10, 51.7, 'Sensitivity', color = 'k', rotation=-28, fontsize=18)
        plt.grid(lw=0.0)
        #plt.ylim(50, 55.5)
        plt.xlim(-11., -6.)
        time_window_str = '{:.2e} s, '.format(self.delta_t) if self.transient else 'Time integrated, '
        custom_labs = [Line2D([0], [0], color = 'k', lw=2., label='This analysis (' + time_window_str + '{:.1f} yr.)'.format(self.data_years))]
        if compare:
            custom_labs.append(Line2D([0], [0], color='grey', lw=2., label=comp_str))
        plt.legend(handles=custom_labs, loc=1)
        plt.ylabel(self.lumi_label, fontsize = 22)
        plt.xlabel(self.density_label, fontsize = 22)
        title = self.lumi_str + ', ' + self.evol_str
        plt.title(title)
        #plt.show()

    def rotated_sensitivity_plot_ts(self, log_ts=False, in_ts=True, ts_vs_p=False, compare=False):
        r'''Two dimensional contour plot to show the sensitivity of the analysis
        in the rotated luminosity-density plane, highlighting the energy requirements
        for the diffuse flux

        Parameters
        ----------
            - compare (bool): include constraints from previous analyses
            - log_ts (bool): contour colors from log or linear ts values
            - in_ts (bool): TS value or binomial-p value
            - ts_vs_p (bool): compare TS and binomial-p value construction
        '''
        if in_ts or ts_vs_p:
            if self.background_median_ts is None:
                self.get_overall_background_ts()
            if self.med_TS is None:
                self.get_med_TS()
        if (not in_ts) or ts_vs_p:
            if self.background_median_p is None:
                self.get_overall_background_p()
            if self.med_p is None:
                self.get_med_p()
        fig, ax = plt.subplots(figsize=(8,5), dpi=200)
        fig.set_facecolor('w')
        X, Y = np.meshgrid(self.densities, self.plot_lumis)
        Y *= X #Scale by the densities
        X = np.log10(X); Y = np.log10(Y)
        if in_ts or ts_vs_p:
            plot_vals = self.med_TS if not log_ts else np.log10(self.med_TS) 
        else:
            plot_vals = self.med_p if not log_ts else np.log10(self.med_p) 
        if in_ts:
            if log_ts:
                levels = np.linspace(np.log10(self.background_median_ts), np.log10(self.background_median_ts) + 2., 11)
            else:
                levels = np.linspace(self.background_median_ts, self.background_median_ts * 10., 11)
        else:
            if log_ts:
                levels = np.linspace(-10., np.log10(self.background_median_p), 11)
            else:
                levels = np.linspace(0.0, self.background_median_p, 11)
        cmap = self.cmap if in_ts else ListedColormap(self.cmap.colors[::-1])
        extend = 'max' if in_ts else 'min'
        cs = ax.contour(X, Y, plot_vals, cmap=cmap, levels=levels, 
                        #vmin=-0.5, 
                        extend=extend)
        csf = ax.contourf(X, Y, plot_vals, cmap=cmap, levels=levels, 
                        #vmin=-0.5, 
                        extend=extend)
        cbar = plt.colorbar(csf) 
        if in_ts or ts_vs_p:
            cbar_lab = r'Median Stacked TS' if not log_ts else r'$\log_{10}($Median Stacked TS$)$'
        else:
            cbar_lab = r'Median binom. p' if not log_ts else r'$\log_{10}($Median binom. p$)$'
        cbar.set_label(cbar_lab, fontsize = 18)
        cbar.ax.tick_params(axis='y', direction='out')
        if in_ts or ts_vs_p:
            cs_ts = ax.contour(X, Y, self.lower_10 - self.background_median_ts, colors=['k'], 
                            levels=[0.0], linewidths=2.)
        if (not in_ts) or ts_vs_p:
            linestyle = 'solid' if ts_vs_p else 'dashed'
            cs_ts = ax.contour(X, Y, self.background_median_p - self.lower_10_p, colors=['k'], 
                            levels=[0.0], linewidths=2., linestyles=linestyle)
        xs = np.logspace(-11., -6., 1000)
        ys_max = self.no_evol_energy_density / xs / self.seconds_per_year if self.transient else self.no_evol_energy_density / xs
        ys_min = self.energy_density / xs / self.seconds_per_year if self.transient else self.energy_density / xs
        plt.fill_between(np.log10(xs), np.log10(ys_min*xs), np.log10(ys_max*xs), 
                color = 'm', alpha = 0.3, lw=0.0, zorder=10)
        if compare:
            comp_rho, comp_en, comp_str = self.compare_other_analyses()
            plt.plot(comp_rho, comp_rho+comp_en, color = 'gray', lw=2.) #damn look at that log property
        plt.text(-10, np.log10(np.max(ys_max*xs)*1.1), 'Diffuse', color = 'm', rotation=0, fontsize=18)
        #plt.text(-10, np.log10(np.min(ys_min*xs)*0.2), 'Sensitivity', color = 'k', rotation=0, fontsize=18)
        plt.grid(lw=0.0)
        plt.xlim(-11., -6.)
        plt.ylim(np.log10(np.min(ys_min*xs)*3e-2), np.log10(np.max(ys_max*xs)*2))
        plt.ylabel(self.scaled_lumi_label, fontsize = 22)
        plt.xlabel(self.density_label, fontsize = 22)
        time_window_str = '{:.2e} s, '.format(self.delta_t) if self.transient else 'Time integrated, '
        custom_labs = [Line2D([0], [0], color = 'k', lw=2., label='This analysis (' + time_window_str + '{:.1f} yr.)'.format(self.data_years))]
        if compare:
            custom_labs.append(Line2D([0], [0], color='grey', lw=2., label=comp_str))
        plt.legend(handles=custom_labs, loc=4)
        title = self.lumi_str + ', ' + self.evol_str
        plt.title(title)
        #plt.show()
    
    def get_labels(self):
        r'''Run during initialization to get the correct units 
        for various plots'''
        self.lumi_label = r'$\log_{10}\Big( \frac{\mathcal{E}}{\mathrm{erg}} \Big)$' if self.transient \
                else r'$\log_{10}\Big( \frac{\mathcal{L}}{\mathrm{erg}\;\mathrm{yr}^{-1}} \Big)$'
        self.density_label = r'$\log_{10}\Big( \frac{\dot{\rho}}{ \mathrm{Mpc}^{-3}\,\mathrm{yr}^{-1}} \Big)$' if self.transient \
                else r'$\log_{10}\Big( \frac{\rho}{ \mathrm{Mpc}^{-3}} \Big)$'
        self.cmap = ListedColormap(sns.light_palette((210, 90, 60), input="husl", n_colors=12))
        self.plot_lumis = self.luminosities / self.time_window_per_year
        self.scaled_lumi_label = r'$\log_{10}\Big( \frac{\mathcal{E}\dot{\rho} }{\mathrm{erg}\;\mathrm{Mpc}^{-3}\,\mathrm{yr}^{-1}} \Big)$' if self.transient \
                else r'$\log_{10}\Big( \frac{\mathcal{L} \rho}{\mathrm{Mpc}^{-3} \mathrm{erg}\;\mathrm{yr}^{-1}} \Big)$'
        self.dens_with_units = r'$\rho$ (Mpc$^{-3}$ yr$^{-1}$)' if self.transient else r'$\rho$ (Mpc$^{-3}$)'
        self.dens_units = r'Mpc$^{-3}$ yr$^{-1}$' if self.transient else r'Mpc$^{-3}$'

    def get_med_TS(self):
        r'''Iterate over the parameter space,
        and extract relevant TS values from distributions
        after trials simulating the parameter space have been run'''
        fmt_path = 'ts_dists_{}year_density_{:.2e}_' + self.evol_lumi_str + \
                        '_manual_lumi_{:.1e}' + self.steady_str + '.npy'
        shape = (self.luminosities.size, self.densities.size)             
        med_TS = np.zeros(shape); lower_10 = np.zeros(shape)
        for ii, lumi in enumerate(self.luminosities):
            for jj, dens in enumerate(self.densities):
                test_en = lumi*dens*self.delta_t if self.transient else lumi*dens
                if test_en > self.energy_density*5.:
                    lower_10[ii, jj] = self.ts_fills[1]
                    med_TS[ii, jj] = self.ts_fills[1]
                elif test_en < self.energy_density*1e-4:
                    lower_10[ii, jj] = self.background_lower_10_ts
                    med_TS[ii, jj] = self.background_median_ts
                else:
                    try:
                        trials = np.load(self.ts_path \
                            + fmt_path.format(self.data_years, dens, lumi))
                        lower_10[ii, jj] = np.percentile(trials[0], 10.)
                        med_TS[ii, jj] = np.median(trials[0])
                    except IOError, e:
                        lower_10[ii, jj] = np.nan
                        med_TS[ii, jj] = np.nan
        med_TS = np.where(np.isnan(med_TS), self.background_median_ts, med_TS)
        lower_10 = np.where(np.isnan(lower_10), self.background_lower_10_ts, lower_10)
        self.med_TS = med_TS
        self.lower_10 = lower_10

    def get_med_p(self):
        r'''Iterate over the parameter space,
        and extract relevant Binomial-p values values from distributions
        after trials simulating the parameter space have been run'''
        fmt_path = 'ts_dists_{}year_density_{:.2e}_' + self.evol_lumi_str + \
                        '_manual_lumi_{:.1e}' + self.steady_str + '.npy'
        shape = (self.luminosities.size, self.densities.size)             
        med_p = np.zeros(shape); lower_10_p = np.zeros(shape)
        for ii, lumi in enumerate(self.luminosities):
            for jj, dens in enumerate(self.densities):
                test_en = lumi*dens*self.delta_t if self.transient else lumi*dens
                if test_en > self.energy_density*5.:
                    lower_10_p[ii, jj] = 1e-10
                    med_p[ii, jj] = 1e-10
                elif test_en < self.energy_density*1e-4:
                    lower_10_p[ii, jj] = self.background_lower_10_p
                    med_p[ii, jj] = self.background_median_p
                else:
                    try:
                        trials = np.load(self.ts_path \
                            + fmt_path.format(self.data_years, dens, lumi))
                        lower_10_p[ii, jj] = np.percentile(trials[2], 90.)
                        med_p[ii, jj] = np.median(trials[2])
                    except IOError, e:
                        lower_10_p[ii, jj] = np.nan
                        med_p[ii, jj] = np.nan
        med_p = np.where(np.isnan(med_p), self.background_median_p, med_p)
        lower_10_p = np.where(np.isnan(lower_10_p), self.background_lower_10_p, lower_10_p)
        self.med_p = med_p
        self.lower_10_p = lower_10_p

    def one_dim_ts_distributions(self, only_gold=False, in_ts = True, log_ts=True):
        r'''Assuming that the diffuse flux is saturated,
        show brazil band plot that scans over density
        and shows the TS distributions'''
        ts_or_p = 0 if in_ts else 1
        ts_inds = (0, 2) if not only_gold else (1, 3)
        levels = []; dens = []
        for density in self.densities:
            try:
                ts = np.load(self.ts_path + \
                    'ts_dists_{}year_density_{:.2e}_'.format(self.data_years, density) + \
                    self.evol_lumi_str + self.steady_str + '.npy')
            except Exception, e:
                continue
            dens.append(density)
            levels.append(np.percentile(ts[ts_inds[ts_or_p]], [5, 25, 50, 75, 95]))
        levels = np.array(levels).T
        fig = plt.figure(dpi=150, figsize=(8,5))
        fig.set_facecolor('w')
        plt.fill_between(dens, levels[0], levels[-1], alpha = 0.5, 
                            color = sns.xkcd_rgb['light navy blue'], 
                            linewidth = 0.0, label = 'Central 90\%')
        plt.fill_between(dens, levels[1], levels[-2], alpha = 0.75, 
                            color = sns.xkcd_rgb['light navy blue'], 
                            linewidth = 0.0, label = 'Central 50\%')
        plt.plot(dens, levels[2], color = sns.xkcd_rgb['light navy blue'])
        plt.title("{}, {}".format(self.lumi_str, self.evol_str))
        plt.xlabel(self.dens_with_units)
        ylab = 'TS' if in_ts else 'Binomial p'
        plt.ylabel(ylab)
        plt.xscale('log')
        loc = 1 if in_ts else 4
        plt.legend(loc=loc, frameon=False)
        if log_ts:
            plt.yscale('log')
        #plt.show()

    def ts_and_ps_plot(self, only_gold=False, log_ts=True):
        r'''Make TS distributions for density, luminosity 
        pairs that saturate the diffuse flux'''
        ts_inds = (0, 2) if not only_gold else (1, 3)
        fig, axs = plt.subplots(ncols=2, nrows=1, dpi=200, sharey=True, figsize=(10,4))
        plt.subplots_adjust(wspace=0.08)
        for density in self.densities[::4]:
            try:
                ts = np.load(self.ts_path + 'ts_dists_{}year_density_{:.2e}_'.format(self.data_years, density)
                                + self.evol_lumi_str + self.steady_str + '.npy')
            except IOError, e:
                continue
            ts_bins = np.logspace(-1., 2., 31) if log_ts else np.linspace(0., 15., 31)
            axs[0].hist(ts[ts_inds[0]], bins = ts_bins, label = r'$\rho = $' 
                    + '{:.1e}'.format(density) + ' ' + self.dens_units, 
                    histtype='step', lw=2.5, weights = [1./len(ts[0])]*len(ts[0]))
            axs[1].hist(ts[ts_inds[1]], bins = np.logspace(-20., 0, 31), label = r'$\rho = $' 
                    + '{:.1e}'.format(density) + ' ' + self.dens_units, 
                    histtype='step', lw=2.5, weights = [1./len(ts[2])]*len(ts[2]))
        fig.suptitle(self.lumi_str + '\n'+ self.evol_str, y=1.03)
        axs[0].set_ylabel('Probability')
        axs[0].set_xlabel('TS')
        axs[1].set_xlabel(r'Binomial $p$')
        if log_ts:
            axs[0].set_xscale('log')
        axs[1].set_xscale('log')
        axs[0].set_yscale('log'); axs[1].set_yscale('log')
        axs[0].set_ylim(8e-3, 7e-1); axs[1].set_ylim(8e-3, 7e-1)
        axs[1].legend(loc=(1.01, 0.1))
        #plt.show()

    def get_overall_background_ts(self, n_trials=5000):
        r'''Sample alert event background distributions
        to get the overall stacked background ts distribution'''
        if self.background_median_ts is not None:
            return self.background_median_ts
        if self.sigs is None:
            self.sigs = self.load_signalness_array()
        bg_trials = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/bg/'
        TSs = []
        for ind in range(len(skymap_files)):
            if self.transient and self.delta_t == 1000.:
                problem_inds = [60, 79, 228]
            elif self.transient:
                problem_inds = [60]
            else:
                problem_inds = [13, 32, 60, 83, 143, 147]
            if ind in problem_inds:
                continue
            else:
                smeared_str = 'smeared/' if self.smeared else 'norm_prob/'
                if self.transient:
                    trials_file = glob(bg_trials + smeared_str 
                                + 'index_{}_*_time_{:.1f}.pkl'.format(ind, self.delta_t))[0]
                    trials = np.load(trials_file)
                    ts = np.random.choice(trials['ts_prior'], size=n_trials)
                else:
                    try:
                        trials_files = glob(bg_trials + smeared_str 
                                    + 'index_{}_steady_seed_*.pkl'.format(ind))
                        trials = []
                        for f in trials_files:
                            trials.extend(np.load(f)['TS'])
                        ts = np.random.choice(np.array(trials), size=n_trials)
                    except:
                        print("NO STEADY TRIALS FOR INDEX {}".format(ind))
                TSs.append(ts)
        TSs = np.array(TSs)
        stacked_ts = np.multiply(TSs, self.sigs[:, np.newaxis])
        stacked_ts = np.sum(stacked_ts, axis=0) / (stacked_ts.shape[0] - 1.) #-1 because we skip one of the maps
        self.background_median_ts = np.median(stacked_ts)
        self.background_lower_10_ts = np.percentile(stacked_ts, 10.)
        self.stacked_ts = stacked_ts
        return self.background_median_ts

    def get_overall_background_p(self, n_trials=5000):
        r'''Sample alert event background distributions
        to get the overall stacked background binomial-p value distribution'''
        bg_trials = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/analysis_trials/bg/'
        if self.sigs is None:
            self.sigs = self.load_signalness_array()
        pvals = []
        for ind in range(len(skymap_files)):
            if self.transient and self.delta_t == 1000.:
                problem_inds = [60, 79, 228]
            elif self.transient:
                problem_inds = [60]
            else:
                problem_inds = [13, 32, 60, 83, 143, 147]
            if ind in problem_inds:
                continue
            else:
                smeared_str = 'smeared/' if self.smeared else 'norm_prob/'
                if self.transient:
                    trials_file = glob(bg_trials + smeared_str 
                                + 'index_{}_*_time_{:.1f}.pkl'.format(ind, self.delta_t))[0]
                    trials = np.load(trials_file)
                    cdf = np.cumsum(np.sort(np.array(trials['ts_prior'])))
                else:
                    trials_files = glob(bg_trials + smeared_str 
                                + 'index_{}_steady_seed_*.pkl'.format(ind))
                    trials = []
                    for f in trials_files:
                        trials.extend(np.load(f)['TS'])
                    cdf = np.cumsum(np.sort(np.array(trials)))
                inds = np.linspace(0., 1., len(cdf))
                inds = np.where(inds==0., np.min(inds[inds != 0.]), inds)[::-1]
                ps = np.where(cdf != 0.0, inds, 1.0)
                pvals.append(np.random.choice(ps, size = n_trials))
        pvals = np.array(pvals)
        background_binomial = []; counter = 0;
        for realization in pvals.T:
            counter += 1
            realization = np.sort(realization)
            obs_p = 1.
            for i, p in enumerate(realization):
                tmp = st.binom_test(i+1, len(realization), p, alternative='greater')
                if tmp < obs_p and tmp != 0.0:
                    if tmp == 0.0:
                        print("WHY DOES THE BINOMIAL VALUE EQUAL ZERO")
                    obs_p = tmp
            background_binomial.append(obs_p)
        background_binomial = np.array(background_binomial)
        binomial_median = np.percentile(background_binomial, 50.)
        binomial_lower_10 = np.percentile(background_binomial, 90.)
        self.background_median_p = np.percentile(background_binomial, 50.)
        self.background_lower_10_p = np.percentile(background_binomial, 90.)
        self.stacked_p = background_binomial
        return self.background_median_p

    def inject_and_fit_dens_lumi_plot(self, dens, lumi, in_ts=True, upper_limit=False):
        r'''Assume a certain density and luminosity, 
        inject it, and see what confidence intervals we 
        can construct
        
        Parameters:
        -----------
            - dens (float): Density of sources
            - lumi (float): Luminosity of sources
        '''
        fmt_path = 'ts_dists_{}year_density_{:.2e}_' + self.evol_lumi_str + \
                        '_manual_lumi_{:.1e}' + self.steady_str + '.npy'
        trials = np.load(self.ts_path + fmt_path.format(self.data_years, dens, lumi))
        trials = trials[0] if in_ts else trials[2]
        unblinded_val = np.random.choice(trials)
        self.inject_and_fit_TS_plot(unblinded_val, in_ts=in_ts, show=False, title=False, upper_limit=upper_limit)
        plt.scatter(np.log10(dens), np.log10(lumi / self.time_window_per_year), 
                        marker='*', color = 'k', s=100)
   

    def inject_and_fit_TS_plot(self, unblinded_val, in_ts=True, show=True, 
                    title=True, upper_limit=False):
        r'''Assume a certain unblinded TS value 
        or binomial p-value and see what confidence intervals we 
        can construct
        
        Parameters:
        -----------
            - unblinded_val (float): Unblinded TS or binomial p-value
            - in_ts (float): Use stacked TS construction instead of binomial p-value
        '''
        fig, ax = plt.subplots(dpi=200)
        X, Y = np.meshgrid(np.log10(self.densities), np.log10(self.plot_lumis))
        cont = self.TS_constraints(unblinded_val, in_ts=in_ts, upper_limit=upper_limit)
        levels = [90., 100.] if upper_limit else [0., 50., 90.]
        csf = ax.contourf(X, Y, cont, cmap=self.cmap, levels = levels)
        xs = np.logspace(-11., -6., 1000)
        ys_max = self.no_evol_energy_density / xs / self.seconds_per_year if self.transient else self.no_evol_energy_density / xs
        ys_min = self.energy_density / xs / self.seconds_per_year if self.transient else self.energy_density / xs
        plt.fill_between(np.log10(xs), np.log10(ys_min), np.log10(ys_max), 
                color = 'm', alpha = 0.3, lw=0.0, zorder=10)
        if not upper_limit:
            legend_elements = [Patch(facecolor=csf.cmap.colors[0], label='50\%'),
                        Patch(facecolor=csf.cmap.colors[-2], label='90\%')]
            ax.legend(handles=legend_elements, loc=3)
        ax.set_ylim(np.log10(ys_min.min()*0.5), np.log10(ys_max.max()*2.))
        ax.grid(lw=0.0)
        ax.set_ylabel(self.lumi_label, fontsize = 22)
        ax.set_xlabel(self.density_label, fontsize = 22)
        if in_ts:
            title_str = 'Observed TS={:.1e}'.format(unblinded_val)
        else:
            title_str = 'Observed binom. p={:.1e}'.format(unblinded_val)
        if title:
            plt.title(title_str)
        if show:
            plt.show()

    def TS_constraints(self, obs_val, in_ts = True, upper_limit=False):
        r'''Based on the observed value, get the 
        frequentist confidence intervals

        Parameters:
        -----------
            - obs_val (float): Observed TS of binomial p-value
            - in_ts (bool): Stacked TS or binomial p-value construction
            - upper_limit (bool): If true, return as value compatible with upper limit
        '''
        fmt_path = 'ts_dists_{}year_density_{:.2e}_' + self.evol_lumi_str + \
                        '_manual_lumi_{:.1e}' + self.steady_str + '.npy'
        ts_p_ind = 0 if in_ts else 2
        containment = np.zeros((self.luminosities.size, self.densities.size)); 
        for ii, lumi in enumerate(self.luminosities):
            for jj, dens in enumerate(self.densities):
                test_en = lumi*dens*self.delta_t if self.transient else lumi*dens
                if test_en > self.energy_density*5.:
                    containment[ii, jj] = 100.
                elif test_en < self.energy_density*1e-4:
                    containment[ii, jj] = 0.
                else:
                    try:
                        trials = np.load(self.ts_path \
                            + fmt_path.format(self.data_years, dens, lumi))
                        containment[ii, jj] = sp.stats.percentileofscore(trials[ts_p_ind], obs_val)
                    except IOError, e:
                        containment[ii, jj] = 0.
        if upper_limit and in_ts:
            containment = (50. - containment)*2.
        elif upper_limit and not in_ts:
            containment = (50. + containment)*2.
        else:
            containment = np.abs(50. - containment)*2.
        return containment

    def compare_other_analyses(self):
        r'''Get the sensitivities / upper limits
        from previous IceCube analyses'''
        if self.transient:
            nora_comparison = {}
            for key in ['GRB_lims', 'GRB_diffuse', 'CCSN_lims', 'CCSN_diffuse']:
                nora_tmp = np.genfromtxt('/data/user/apizzuto/fast_response_skylab/alert_event_followup/effective_areas_alerts/Nora_{}.csv'.format(key),
                                        delimiter=', ')
                nora_comparison[key] = nora_tmp
            for key in ['GRB_lims']:
                tmp = np.log10(zip(*nora_comparison[key]))
                #plt.plot(tmp[0], tmp[1], color = 'grey')
            return tmp[0], tmp[1], 'Multiplets (100 s, 5 yr.)'
        else:
            ps_pap = np.genfromtxt('/data/user/apizzuto/fast_response_skylab/alert_event_followup/effective_areas_alerts/point_source_paper_lims.csv',
                            delimiter=', ')
            tmp = np.array(zip(*ps_pap))
            lums = tmp[0] * self.no_evol_energy_density/self.energy_density #testing diff. spectrum, this is an approximation rn
            return np.log10(tmp[1]), np.log10(lums), '8 yr. point source'

    def load_signalness_array(self):
        r''' Load the signalness by index list, apply
        appropriate masks depending on which analysis is being run'''
        sigs_all = np.load('/data/user/apizzuto/fast_response_skylab/alert_event_followup/effective_areas_alerts/sigs_by_ind.npy')
        if self.transient and self.delta_t == 1000.:
            msk_inds = np.array([60, 79, 228])
        elif self.transient:
            msk_inds = np.array([60])
        else:
            msk_inds = np.array([13, 32, 60, 83, 143, 147])
        sigs = np.delete(sigs_all, msk_inds)
        return sigs

    def upper_limit_plot(self):
        pass

    def fit_coverage_plot(self, dens, lumi):
        pass

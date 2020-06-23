import numpy as np
from glob import glob
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
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

energy_density = {'transient': {'HB2006SFR': 4.8038e51, 
                        'MD2014SFR': 6.196e51,
                        'NoEvolution': 1.8364e52},
                'steady': {'HB2006SFR': 7.6735e43, 
                        'MD2014SFR': 9.897e+43,
                        'NoEvolution': 2.93335e44}}

class UniversePlotter():
    r'''
    Tool to make some helpful plots
    '''
    def __init__(self, delta_t, data_years, lumi, evol, **kwargs):
        self.delta_t = delta_t
        self.transient = True if self.delta_t is not None else False
        self.data_years = data_years
        self.lumi = lumi
        self.evol = evol
        self.lumi_str = {'SC': 'Standard Candle', 'LG': 'LogNormal'}[self.lumi]
        self.evol_str = {'HB2006SFR': 'Hopkins and Beacom 2006 SFR',
                            'MD2014SFR': 'Madau and Dickinson 2014 CSFH'}[self.evol]
        self.ts_path = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/ts_distributions/'
        self.steady_str = '' if self.transient else '_steady'
        key = 'transient' if self.transient else 'steady'
        self.energy_density = energy_density[key][evol]
        self.no_evol_energy_density = energy_density[key]['NoEvolution']
        self.evol_lumi_str = 'evol_{}_lumi_{}{}'.format(self.evol, self.lumi, self.steady_str)
        self.densities = np.logspace(-11., -6., 21)
        self.luminosities = None #NEED SOME LOGIC HERE

    def initialize_plotter(self):
        '''SIMULATE A TRANSIENT UNIVERSE JUST TO GET THE LUMINOSITY'''
        pass

    def two_dim_sensitivity_plot(self):
        pass

    def one_dim_ts_distributions(self):
        pass

    def ts_and_ps_plot(self):
        for density in np.logspace(-11., -6., 11)[:]:
            try:
                ts = np.load(self.ts_path + 'tsdists_{}year_density_{:.2e}'.format(self.data_years, dens)
                                + self.evol_lumi_str + '.npy')
            except:
                continue
        pass

    def inject_and_fit(self, dens, lumi):
        pass

    def upper_limit(self, TS):
        pass

    def compare_other_analyses(self):
        pass

    def fit_coverage_plot(self, dens, lumi):
        pass
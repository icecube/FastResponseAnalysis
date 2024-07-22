r''' General Fast Response Analysis Class.

    Author: Alex Pizzuto
    Date: 2021

    One of two variants to include mutiple datasets.
    Adding a new class that holds multiple FastResponseAnalysis as a container,
    adds the appropriate LLH, injector, and adapted methods.

    '''

from abc import abstractmethod
import os, sys, time, subprocess
import pickle, dateutil.parser, logging, warnings
from argparse import Namespace
from copy import deepcopy

import h5py
import healpy                 as hp
import numpy                  as np
import seaborn                as sns
import matplotlib             as mpl
import matplotlib.pyplot      as plt
import numpy.lib.recfunctions as rf
from astropy.time           import Time
from scipy.special          import erfinv
from matplotlib.lines       import Line2D

from skylab.datasets        import Datasets
from skylab.llh_models      import EnergyLLH
from skylab.priors          import SpatialPrior
from skylab.ps_injector     import PointSourceInjector
from skylab.ps_llh          import PointSourceLLH, MultiPointSourceLLH
from skylab.ps_injector     import PriorInjector
from skylab.spectral_models import PowerLaw 
from skylab.temporal_models import BoxProfile, TemporalModel
# import meander

from . import web_utils
from . import sensitivity_utils
from . import plotting_utils
from .reports import FastResponseReport
from .FastResponseAnalysis import FastResponseAnalysis, PointSourceFollowup

mpl.use('agg')
current_palette = sns.color_palette('colorblind', 10)
logging.getLogger().setLevel(logging.ERROR)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

class MultiFastResponseAnalysis(FastResponseAnalysis):
    """
    Instead of sprinkling if-else statements around, could leave
    the LLH contained within the original FastResponseAnalysis.
    Analogous to MultiPointSourceLLH, this one then is a container
    for one FRA per dataset, plus a MultiPointSourceLLH, and the injector.
    Will involve some duplicated code, unless the methods are split up more,
    or both FRA and MultiFRA get a new base class (analogous to BaseLLH). 
    But it separates code more cleanly, and can override inherited methods
    instead of having two implementations with an if-else statement inbetween.
    Requires some extra arguments to FastResponseAnalysis methods.
    """

    # inherits some default config attributes
    def __init__(self,  *args, **kwargs):
        """The arguments dataset and index override the defaults 
        that may have been set by child classes.

        Same args and kwargs as FastResponseAnalysis
        
        """

        # basic config of this instance, and initialize_llh
        super().__init__(*args, **kwargs)

        # save spectrum and time profile to help broadcast
        self.spectrum = None
        if self._fix_index:
            self.spectrum = PowerLaw(A=1, gamma=self.index, E0=1000.)
        self.time_profile = BoxProfile(self.start, self.stop)

    @abstractmethod
    def initialize_analyses(self, *args, **kwargs):
        pass

    def initialize_llh(self, skipped=None, scramble=False):
        if self._verbose:
            print("Initializing MultiPointSourceLLH in Skylab")
        
        # kwargs for the BaseLLH instance
        # TODO maybe this should be an attribute,
        # so config can be changed in one place for both FRA snd MultiFRA?
        base_kwargs = dict(
            nsource=1.,                    # seed for nsignal fit
            nsource_bounds=(0., 1e3),      # bounds on fitted ns
            ncpu=self._ncpu,               # use 10 CPUs when computing trials
            seed=self.llh_seed,
        )

        multi_llh = MultiPointSourceLLH(**base_kwargs)
        for enum, fra in enumerate(self.analyses):
            fra.llh.do_trials_seed = 1
            multi_llh.add_sample(fra._dataset, fra.llh)
        return multi_llh
        
    def remove_event(self, exp, dset, skipped):
        # should this ever be called from this instance?
        # need to save event in self.save_items
        raise NotImplemented('')
    
    # These properties are initialized for the base class
    # could do something meaningful
    @property
    def skipped_event(self):
        return self._skipped_event
    # TODO could be a property
    # that retrieves unique skipped_event from constituent analyses
    @skipped_event.setter
    def skipped_event(self, x):
        self._skipped_event = x
    # TODO make something meaningful of this

    

    @property
    def exp(self):
        '''Dictionary or flat array like llh_exp? used in remove_event'''
        raise NotImplemented('')

    @exp.setter
    def exp(self, x):
        self._exp = x

    @property
    def mc(self): # and livetime, grl, sinDec_bins, energy_bins... only used in LLH
        '''Dictionary? i.e. self.llh.mc'''
        raise NotImplemented('')

    @property
    def livetime(self):
        '''Dictionary?'''
        raise NotImplemented('')

    @property
    def dset(self):
        '''Dictionary? used in remove_event'''
        raise NotImplemented('')


    @property
    def datasets(self):
        """Returns the datasets used"""
        return [_a._dataset for _a in self.analyses]
    
    @property
    def analyses(self):
        '''Returns the constituent followup instances'''
        return self._analyses
    @analyses.setter
    def analyses(self, x):
        self._analyses = x
        # Feature to come for easier comparison: 
        # re-initialize LLH from existing FRA's
        # if isinstance(x, list):
        #     self._analyses = x
        # elif isinstance(x, FastResponseAnalysis):
        #     self._analyses = [x]
        # else:
        #     raise TypeError(f'trying to set analyses with {type(x)}')

    @property
    def llh_exp(self):
        """Returns a flat array of experimental data loaded across the used dataset(s) loaded into the LLH,
        limited to fields used for plotting and amended with a field for the dataset."""
        if not hasattr(self, '_llh_exp'):
            exp = {enum:_llh.exp for enum, _llh in self.llh._samples.items()}
            merged_dtype = np.dtype([('run', '<i4'), ('event', '<i4'), ('time', '<f4'), ('ra', '<f4'), ('dec', '<f4'), ('sigma', '<f4'), ('enum', '<i4')])
            for enum, _exp in exp.items():
                _exp = rf.drop_fields(_exp, [_field for _field in _exp.dtype.names if _field not in merged_dtype.names])
                _exp = rf.append_fields(_exp, 'enum', np.full( _exp.size, enum))
                _exp = _exp.astype(merged_dtype)
                exp[enum] = _exp
            self._llh_exp = np.concatenate([exp[enum] for enum in exp])
        return self._llh_exp
    
    # TODO: smarter way to change label
    def plot_skymap(self, **kwargs):
        return super().plot_skymap(label='Event', **kwargs)

class MultiPointSourceFollowup(PointSourceFollowup, MultiFastResponseAnalysis):

    def __init__(self, *args, **kwargs):

        # first, construct constituent analyses with some arguments
        self.initialize_analyses(*args, **kwargs)    
        super().__init__(*args, **kwargs)

    def __str__(self):
        string = super().__str__().rstrip()
        string = string.rstrip()
        string += 'Datasets:\n'
        string += '\n'.join(self.datasets)
        string += '\n\n'
        return string

    def initialize_analyses(self, *args, **kwargs):
        # initialize individual dataset LLH
        # needs to be here as they will have the same constructor signature
        self.analyses = []
        for _followup in self._followups:
            _kwargs = deepcopy(kwargs)
            _kwargs['save'] = False # this single FRA does not save output
            # other class attributes one may want to broadcast
            for attr in ['_verbose',
                         '_index',
                         '_float_index',
                         '_fix_index',
                         '_index_range',
                         '_llh_seed',
                         '_nb_days',
                         '_ncpu',
                         ]:
                setattr(_followup, attr, getattr(self, attr))
            # initialize with analysis specific args and kwargs
            # such as source properties, or override settings
            print(args, kwargs)
            # get tstart, tstop, ra, dec, and other config
            _analysis = _followup(*args, **_kwargs)
            self.analyses.append(_analysis)

    def initialize_injector(self, e_range=(0., np.inf)):
        inj = PointSourceInjector(
            gamma = self.index, 
            E0 = 1000., 
            e_range=e_range)
        temporal_model = {enum:_llh.temporal_model for enum,_llh in self.llh._samples.items()}
        inj.fill(
            self.dec,
            self.llh.exp,
            self.llh.mc,
            self.llh.livetime,
            temporal_model=temporal_model)
        self.inj = inj
        self.save_items['E0'] = self.inj.E0

    def find_coincident_events(self):
        '''Retrieve events with spatial x energy x temporal weight>10
        from all component analyses, store in self.coincident_events
        as a list of dictionaries with columns for the report tables,
        ['run', 'event', 'ra', 'dec', 'sigma', 'logE', 'time']
        added Delta Psi, spatial weight, energy weight, sample enum.
        '''
        # TODO ideally samples should not overlap, however GFU and Greco do. 
        # decide whether to require de-duplication at analysis level,
        # or in this table (and record tuples in enum column).
        for enum, _ana in enumerate(self.analyses):
            _ana.find_coincident_events(ns_params = self.ns_params)
            for _event in _ana.coincident_events:
                _event['enum'] = enum
                _event['dataset'] = _ana.dataset
        self.coincident_events = sum([_ana.coincident_events for _ana in self.analyses], [])
        self.save_items['coincident_events'] = self.coincident_events

    def make_dNdE(self):
        # easier if calculation and plotting were separate methods...
        # but how to do it is already laid out in FRA_ifelse
        r"""Make an E^-2 or E^-2.5 dNdE with the central 90% 
        for the most relevant declination band 
        (+/- 5 deg around source dec)
        for multiple datasets
        """
        low5 = []
        high5 = []
        fig, ax = plt.subplots(figsize = (8,5))
        fig.set_facecolor('white')
        for enum in self.llh._samples:
            llh = self.llh._samples[enum]
            dataset = self.datasets[enum].replace('_', ' ')

            dec_mask_1 = llh.mc['dec'] > self.dec - (5. * np.pi / 180.)
            dec_mask_2 = llh.mc['dec'] < self.dec + (5. * np.pi / 180.)
            dec_mask_3, dec_mask_4 = None, None
            dec_mask = dec_mask_1 * dec_mask_2
            
            
            lab = dataset
            color = f'C{enum+1}'
            delta_gamma = -1. * self.index + 1.
            a = plt.hist(llh.mc['trueE'][dec_mask], bins = np.logspace(1., 8., 50), 
                    weights = llh.mc['ow'][dec_mask] * np.power(llh.mc['trueE'][dec_mask], delta_gamma) / llh.mc['trueE'][dec_mask], 
                    histtype = 'step', linewidth = 2., color = color, label = lab)
            cdf = np.cumsum(a[0]) / np.sum(a[0])
            low_5 = np.interp(0.05, cdf, a[1][:-1])
            median = np.interp(0.5, cdf, a[1][:-1])
            high_5 = np.interp(0.95, cdf, a[1][:-1])
        
            plt.axvspan(low_5, high_5, color = color, alpha = 0.25, label="Central 90\%")
            lab = 'Median'
            plt.axvline(median, c = color, alpha = 0.75, label = lab)
            low5.append(low_5)
            high5.append(high_5)
            print('{} 90% Central Energy Range: {}, {} GeV'.format(dataset, round(low_5), round(high_5)))
            
        plt.yscale('log')
        plt.xscale('log')
        plt.grid(which = 'major', alpha = 0.25)
        plt.xlabel('Energy (GeV)', fontsize = 24)

        plt.xlim(1e1, 1e8)
        plt.legend(loc=4, fontsize=18)
        plt.savefig(self.analysispath + '/central_90_dNdE.png',bbox_inches='tight')

        self.low5 = low5
        self.high5 = high5
        self.save_items['energy_range'] = tuple(zip(self.low5, self.high5))
        return

    def write_circular(self):
        raise NotImplemented('This method is not implemented for the parent class, either.')




    

    
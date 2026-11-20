r''' Multisample Fast Response Analysis Class.

    Author: Christoph Raab
    Date: 2024

    Adding a new class that holds multiple FastResponseAnalysis as a container,
    adds the appropriate LLH, injector, and adapted methods.
    Since FastResponseAnalysis has methods that are (made) dataset-agnostic,
    it can inherit these as well.
    '''

from abc import abstractmethod
import os, sys, time, subprocess
import pickle, dateutil.parser, logging, warnings
from argparse import Namespace
from copy import deepcopy
from collections import defaultdict

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

from . import web_utils
from . import sensitivity_utils
from . import plotting_utils
from .reports import FastResponseReport
from .FastResponseAnalysis import FastResponseAnalysis, PriorFollowup, PointSourceFollowup
from .GWFollowup import GWFollowup
from .AlertFollowup import AlertFollowup, CascadeFollowup, TrackFollowup

mpl.use('agg')
current_palette = sns.color_palette('colorblind', 10)
logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

class MultiFastResponseAnalysis(FastResponseAnalysis):
    """
    This base class supports multi-sample fast response analyses.
    Each sample's LLH is contained within the original FastResponseAnalysis,
    or one of its child classes. This class then
    - contains the classes in self._followups,
    - constructs the instances in self.analyses,
    - combines their data into self.llh_exp,
    - combines their likelihoods into a MultiPointSourceLLH.
    This classes children (MultiPriorFollowup, MultiPointSourceFollowup) then
    - construct the injector
    - override plotting methods where needed,
    (while some are already made to work with single- or multi-sample follow-ups).
    """

    # attributes that will be set for a specific followup configuration inheriting from this base class
    _followups = []

    # attributes that are identical for each analysis produced by one followup class
    # will be broadcast to the constituent analyses
    static_attributes = []

    # inherits some default config attributes
    def __init__(self,  *args, **kwargs):
        """Called from the children
        - PriorFollowup.__init__
        - PointSourceFolloup.__init__

        Same args and kwargs as FastResponseAnalysis, and sets the same attributes,
        only self.llh will be a MultiPointSourceLLH from MultiFastResponseAnalysis.initialize_llh.

        Parameters
        ----------
        name : str
            name of the source
        tstart : float
            start of the follow-up window, [MJD]
        tstop : float
            stop of the follow-up window [MJD]
        skipped : list, optional
            Skipped events, by default None
        seed : int, optional
            LLH RNG seed, by default 1
        outdir : str, optional
            Output base directory, by default will use ${FAST_RESPONSE_OUTPUT} or CWD
        save : bool, optional
            Save analysis results, by default True
        extension : float, optional
            Point source extension, by default None
        index : float, optional
            Spectral index to override class attribute self._index of the child class
        fix_index : bool, optional
            Keep spectrum fixed, to override self._fix_index of the child class
        dataset : str, optional
            Dataset descriptor attribute, to override self._dataset of the child class

        To help broadcast the hypothesis to the constituents it also sets
        self.spectrum (if fixed power law spectrum)
        self.time_profile (common box profile)
        """
        logger.debug('MultiFastResponseAnalysis.__init__')

        # basic config of this instance, and initialize_llh
        super().__init__(*args, **kwargs)

        # save spectrum and time profile to help broadcast
        self.spectrum = None
        if self._fix_index:
            self.spectrum = PowerLaw(
                A=1, # normalization factor for dN/dE (flux) at E0 [GeV^-1cm^-2s^-1]
                gamma=self.index, # spectral index
                E0=1000., # pivot energy [GeV]
                )
        self.time_profile = BoxProfile(self.start, self.stop)

    def initialize_analyses(self, *args, **kwargs):
        """Constructs the constituent followups, turning
        self._followups (classes) => self.analyses (instances)
        with matching settings broadcast via *args and **kwargs.
        """
        # initialize individual dataset LLH
        # needs to be here as they will have the same constructor signature
        self.analyses = []
        for _followup in self._followups:
            #if not issubclass(type(_followup), type(self)):
            #    raise TypeError(f'Trying to construct a {type(self)} with a {type(_followup)}')
            _kwargs = deepcopy(kwargs)
            _kwargs['save'] = False # this single FRA does not save output
            # other class attributes one may want to broadcast
            for attr in self.static_attributes:
                setattr(_followup, attr, getattr(self, attr))
            # initialize with analysis specific args and kwargs
            # such as source properties, or override settings
            
            # get tstart, tstop, ra, dec / skymap, and other config
            _analysis = _followup(*args, **_kwargs)
            self.analyses.append(_analysis)

    def initialize_llh(self, skipped=None, scramble=False):
        """Once self.analyses are constructed, this constructs
        self.llh = the MultiPointSourceLLH combining their LLH's.
        """
        if self._verbose:
            print("Initializing MultiPointSourceLLH in Skylab")
        
        # kwargs for the BaseLLH instance
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
        raise NotImplementedError('remove_event not appropriate for MultiFastResponseAnalysis')
    
    @property
    def mc(self):
        """Dictionary {dataset enum : dataset mc array}"""
        return self.llh.mc

    @property
    def livetime(self):
        """Dictionary {dataset enum : dataset livetime}"""
        return self.llh.livetime

    @property
    def datasets(self):
        """List [used dataset names]"""
        return [_a._dataset for _a in self.analyses]
    
    @property
    def analyses(self):
        '''List [followup instances]'''
        return self._analyses
    @analyses.setter
    def analyses(self, x):
        self._analyses = x
        
    # We share an LLH seed
    @property
    def llh_seed(self):
        """RNG seed, MultiPSLLH propagates it to the PSLLHs"""
        return self._llh_seed

    # when setting a new one, broadcast it
    @llh_seed.setter
    def llh_seed(self, value):
        self._llh_seed = value
        if hasattr(self, "llh"):
            logger.warning("Propagating new RNG seed to LLHs")
            self.llh.set_rng_seed(value)
        

    @property
    def llh_exp(self):
        """Flat array of experimental data loaded across the used dataset(s) loaded into the LLH,
        limited to fields used for plotting and amended with a field for the dataset."""
        if not hasattr(self, '_llh_exp'):
            exp = {enum:_llh.exp for enum, _llh in self.llh._samples.items()}
            for enum, _exp in exp.items():
                exp[enum] = self.unify_exp_array(_exp, enum=enum)
            self._llh_exp = np.concatenate([exp[enum] for enum in exp])
        return self._llh_exp
    
    def plot_skymap(self, **kwargs):
        """Plot skymap of all dataset events with short dataset legend.

        **kwargs as in FastResponseAnalysis.plot_skymap, except labels
        with_contour: bool
            plots the 90% containment contour of a skymap (default False)
        contour_files: string
            text file containing skymap contours to be plotted (default None)
        label_events: bool
            adds a number label to events on skymap (default False)
        distinct_colorbars : bool
            Use distinct color bars for each dataset's time axis (default False)
        """
        labels = []
        for _ds in self.datasets:
            _base = _ds.split('_')[0] # convention: version after underscore
            _base = _base.replace('Online', '') # not necessary for legend
            _base = _base.replace('Greco', 'GRECO') # some prefer this
            labels.append(f'{_base} Event')
        return super().plot_skymap(labels=labels, **kwargs)

    @property
    def dataset_string(self):
        """String listing used datasets for the __repr__"""
        string = ''
        string += 'Datasets:\n'
        string += '\n'.join(self.datasets)
        string += '\n\n'
        return string

class MultiPriorFollowup(PriorFollowup, MultiFastResponseAnalysis):
    """"
    Base class for multi-dataset followups using a spatial prior.
    It broadcasts properties to constituent PriorFollowup's configured as
    a list of classes in the _followup class attribute, then
    combines their likelihoods into a MultiPointSourceLLH which can
    still use the same skylab.base_llh.BaseLLH.scan() method.
    """
    # attributes that are identical for each analysis produced by one followup class
    # will be broadcast to the constituent analyses
    static_attributes = ['_verbose',
                         '_index',
                         '_float_index',
                         '_fix_index',
                         '_index_range',
                         '_llh_seed',
                         '_ncpu',
                         '_pixel_scan_nsigma',
                         '_allow_neg',
                         '_containment',
                         '_nside',
                         ]

    def __init__(self, *args, **kwargs):
        """Construct the MultiPointSourceFollowup,
        1) the constituent followups self.analyses,
        2) its own likelihood self.llh
        3) its precomputed data directories self._bg_dir, self._sens_dir

        Parameters
        ----------
        *args, **kwargs
            Passed to self.initialize_analyses and MultiFastResponseAnalsis.__init__
        """
        logger.debug('MultiPriorFollowup.__init__')
        # first, construct constituent analyses with same arguments
        self.initialize_analyses(*args, **kwargs)

        # then prepare LLH and store prior-specific attributes
        super().__init__(*args, **kwargs)
        # precomputed sensitivity and trial can have a generic path
        # but actually specified by the configured datasets
        self._bg_dir = os.path.join(self._bg_dir, "+".join(self.datasets))
        self._sens_dir = os.path.join(self._sens_dir, "+".join(self.datasets))

    def __str__(self):
        string = super().__str__().rstrip()
        string += '\n'
        string += self.dataset_string
        return string
    
    def initialize_injector(self, e_range=(0., np.inf)):
        """Construct the spatial prior based on the skymap and then
        the multi-dataset PriorInjector, using the dictionaries
        {dataset enum: PointSourceLLH attribute} dictionaries
        already present as MultiPointSourceLLH attributes exp, mc and livetime,
        plus a similar one constructed here for the temporal model.
        Assumes power law with self.index and pivot energy 1 TeV.

        Parameters
        ----------
        e_range : tuple, optional
            True energy range to inject [GeV], by default (0., np.inf)
            Used for differential sensitivities.
        """
        logger.debug("Initializing Prior Injector")
        spatial_prior = SpatialPrior(self.skymap, containment = self._containment, allow_neg=self._allow_neg)
        self.spatial_prior = spatial_prior
        inj = PriorInjector(
            spatial_prior, 
            gamma=self.index, 
            e_range = e_range, 
            E0=1000., # pivot [GeV]
            seed = self.llh_seed)
        temporal_model = {enum:_llh.temporal_model for enum,_llh in self.llh._samples.items()}
        inj.fill(
            self.llh.exp,
            self.llh.mc,
            self.llh.livetime,
            temporal_model=temporal_model)
        self.inj = inj
        self.save_items['E0'] = self.inj.E0
    
    def make_dNdE(self):
        r"""Make an E^-2 or E^-2.5 dNdE with the central 90% 
        for the minimum and maximum declinations on the skymap
        for multiple datasets
        """
        min_dec, max_dec = self.dec_skymap_range()
        low5 = []
        high5 = []
        fig, ax = plt.subplots(figsize = (8,5))
        fig.set_facecolor('white')
        
        # iterate over dataset
        for enum in self.llh._samples:
            llh = self.llh._samples[enum]
            dataset = self.datasets[enum].replace('_', ' ')
            style = plotting_utils.skymap_style[enum]

            energy_range = defaultdict(list)

            # then iterate over the min- and max- declination of the skymap
            for (dec_label, dec, color) in [
                ("min. dec", min_dec, sns.xkcd_rgb['windows blue']),
                ("max. dec", max_dec, sns.xkcd_rgb['dark navy blue']),

            ]:
                label = ""
                if len(self.analyses)==1:
                    label = dec_label
                else:
                    label = " ".join((dec_label, dataset))

                # obtain the quantities for this dataset
                energy_band = plotting_utils.get_energy_band(llh.mc, self.index, dec,
                                                            half_width=5., coverage=0.9)
                # then make the plot, differentiating by linestyle
                plotting_utils.plot_energy_band(**energy_band, color=color, linestyle=style["linestyle"],
                                                label_prefix=label)
                for key in ["low", "high", "median"]:
                    energy_range[key].append(energy_band[key])
            
            low5.append(np.min(energy_range["low"]))
            high5.append(np.max(energy_range["high"]))
            
            
        plt.yscale('log')
        plt.xscale('log')
        plt.grid(which = 'major', alpha = 0.25)
        plt.xlabel('Energy (GeV)', fontsize = 24)

        plt.xlim(1e1, 1e8)
        plt.legend(loc=4, fontsize=18)
        plt.savefig(self.analysispath + '/central_90_dNdE.png',bbox_inches='tight')

        self.low5 = low5 # one entry per sample
        self.high5 = high5 # one entry per sample
        self.energy_range =  tuple(zip(self.low5, self.high5))
        self.save_items['energy_range'] = self.energy_range

class MultiPointSourceFollowup(PointSourceFollowup, MultiFastResponseAnalysis):
    """"
    Base class for multi-dataset followups of a point source.
    It broadcasts properties to constituent PointSourceFollowup's configured
    as a list of classes in the _followup class attribute, then
    combines their likelihoods into a MultiPointSourceLLH which can
    still use the same skylab.base_llh.BaseLLH.fit_source() method.
    """
    # attributes that are identical for each analysis produced by one followup class
    # will be broadcast to the constituent analyses
    static_attributes = ['_verbose',
                         '_index',
                         '_float_index',
                         '_fix_index',
                         '_index_range',
                         '_llh_seed',
                         '_ncpu',
                         ]

    def __init__(self, *args, followups=None, **kwargs):
        """Construct the MultiPointSourceFollowup,
        1) the constituent followups
        2) its own likelihood.

        Parameters
        ----------
        followups : list, optional
            Override self._followups

        *args, **kwargs
            Passed to self.initialize_analyses and MultiFastResponseAnalsis.__init__
        """
        logger.debug('MultiPointSourceFollowup.__init__')
        if followups is not None:
            self._followups = followups

        # first, construct constituent analyses with same arguments
        self.initialize_analyses(*args, **kwargs)

        # then construct LLH
        super().__init__(*args, **kwargs)

    def __str__(self):
        string = super().__str__().rstrip()
        string += self.dataset_string
        return string

    def initialize_injector(self, e_range=(0., np.inf)):
        """Construct the point source injector, which SkyLab makes
        rather transparent since it will use the dictionaries
        {dataset enum: PointSourceLLH attribute} dictionaries
        already present as MultiPointSourceLLH attributes exp, mc and livetime.
        A similar one is constructed here for the temporal model.
        Assumes a power law with self.index and pivot energy 1 TeV.

        Parameters
        ----------
        e_range : tuple, optional
            True energy range to inject [GeV], by default (0., np.inf)
            Used for differential sensitivities.
        """
        inj = PointSourceInjector(
            gamma = self.index, 
            E0 = 1000., # pivot [GeV]
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
        # NOTE ideally samples should not overlap and this is prevented by configuring the dataset
        # if not, duplicates will appear appear here twice
        # (which is a warning that the combined LLH is iffy)

        for enum, _ana in enumerate(self.analyses):
            _ana.find_coincident_events(ns_params = self.ns_params)
            for _event in _ana.coincident_events:
                _event['enum'] = enum
                _event['dataset'] = _ana.dataset
        self.coincident_events = sum([_ana.coincident_events for _ana in self.analyses], [])
        self.save_items['coincident_events'] = self.coincident_events

    def make_dNdE(self):
        r"""Make an E^-2 or E^-2.5 dNdE with the central 90% 
        for the most relevant declination band 
        (+/- 5 deg around source dec)
        for multiple datasets
        """
        low5 = []
        high5 = []
        fig, ax = plt.subplots(figsize = (8,5))
        fig.set_facecolor('white')

        # iterate over dataset
        for enum in self.llh._samples:
            llh = self.llh._samples[enum]
            dataset = self.datasets[enum].replace('_', ' ')
            style = plotting_utils.skymap_style[enum]
            color = sns.xkcd_rgb['windows blue']

            # obtain the quantities for this dataset
            energy_band = plotting_utils.get_energy_band(llh.mc, self.index, self.dec,
                                                         half_width=5., coverage=0.9)
            # then make the plot, differentiating by linestyle
            plotting_utils.plot_energy_band(**energy_band, color=color, linestyle=style["linestyle"],
                                            label_prefix=dataset)
            
            low5.append(energy_band["low"])
            high5.append(energy_band["high"])
            
            
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

    def write_circular(self):
        """Write a circular, which needs to be implemented for specific approved analyses in child classes."""
        raise NotImplementedError('This method is not implemented for the parent class, either.')

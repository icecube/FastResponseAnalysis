r''' General Fast Response Analysis Class.
    Heavily influenced by Raamis Hussain's 
    Realtime GW analysis and Kevin Meagher's 
    original Fast Response analysis.

    Author: Alex Pizzuto
    Date: April, 2019
    '''

import os, sys, time, subprocess, pickle, shutil, dateutil.parser, logging, utils
import healpy as hp
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from astropy.io import fits
from astropy.time import Time
from matplotlib import cm
import scipy as sp
from scipy.optimize       import curve_fit
from scipy.stats          import chi2
from scipy.special import erfinv
from make_ontime_plots import make_rate_plots
from ReportGenerator import ReportGenerator
from skylab.datasets import Datasets
from skylab.llh_models import EnergyLLH
from skylab.priors import SpatialPrior
from skylab.ps_injector import PointSourceInjector
from skylab.ps_llh import PointSourceLLH
from skylab.spectral_models import PowerLaw
from skylab.temporal_models import BoxProfile, TemporalModel
# import meander
############################# Plotting Parameters #############################
mpl.use('agg')
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['mathtext.rm'] = 'Times New Roman'
mpl.rcParams['mathtext.it'] = 'Times New Roman:italic'
mpl.rcParams['mathtext.bf'] = 'Times New Roman:bold'

mpl.rc('font', family='serif', size=12)
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['ytick.major.size'] = 5

current_palette = sns.color_palette('colorblind', 10)
############################# Plotting Parameters #############################

logging.getLogger("skylab.ps_llh.PointSourceLLH").setLevel(logging.ERROR)
logging.getLogger("skylab.ps_injector.PointSourceInjector").setLevel(logging.ERROR)

class FastResponseAnalysis(object):
    r''' Object to do realtime followup analyses of 
        astrophysical transients with arbitrary event
        localization
    '''

    def __init__(self,location,tstart,tstop,**kwargs):
        r'''Constructor

        Parameters:
        -----------
        location: str
            Either the location of a point source as RA, Dec 
            in degrees, or a skymap (passed as a url to a .fits
            file or as a healpy file already on cobalt)
        tstart: str
            Beginning of analysis time window in UTC. 
            Must have the form:
            '2019-03-10 12:00:00' 
        tstop: str
            Beginning of analysis time window in UTC. 
            Must have the form:
            '2019-03-10 12:00:00' 

        '''
        self.name = kwargs.pop("Name", "FastResponseAnalysis")
        
        if location.strip('() ').find(' ') > -1:        
            #Point source passed with ra, dec
            location = location.strip('() ').split(' ')
            self.ra, self.dec = np.deg2rad(float(location[0].rstrip(','))), np.deg2rad(float(location[1]))
            if (self.ra < 0 or self.ra > 2*np.pi) or (self.dec < -np.pi/2. or self.dec > np.pi/2.):
                print("Right ascension and declination are not valid")
                sys.exit()
            self.skymap, self.nside, self.ipix_90 = None, None, None
            self.extension = kwargs.pop("Extension", None)
            if self.extension is not None:
                self.extension = np.deg2rad(float(self.extension))
            self.source_type = "PS"
        else: 
            #Skymap passed as a healpy file  
            try: 
                location = location.strip('() ')
                self.skymap_url = location
                if "steinrob" in location:
                    self.skymap = fits.open(location)[0].data
                elif 'data/ana/' in location:
                    skmap_files = glob('/data/ana/realtime/alert_catalog_v2/2yr_prelim/fits_files/Run13*.fits.gz')
                    skymap_f = fits.open(skymap_fies[alert_ind])
                    skymap_fits = skymap_f[1].data
                    skymap_header = skymap_f[1].header
                    run_id, ev_id = skymap_header['RUNID'], skymap_header['EVENTID']
                    ev_mjd = skymap_header['EVENTMJD']
                    ev_iso = skymap_header['START']
                    signalness = skymap_header['SIGNAL']
                    ev_en = skymap_header['ENERGY']
                    ev_ra, ev_dec = np.radians(skymap_header['RA']), np.radians(skymap_header['DEC'])
                    ev_stream = skymap_header['I3TYPE']
                    self.skymap = np.exp(skymap_fits / 2.) 
                    self.skymap /= self.skymap.sum()
                else:
                    self.skymap = hp.read_map(location, verbose=False)
                if hp.pixelfunc.get_nside(self.skymap)!=256:
                    self.skymap = hp.pixelfunc.ud_grade(self.skymap,256)
                    self.skymap = self.skymap/self.skymap.sum()
                self.nside = hp.pixelfunc.get_nside(self.skymap)
                self.ipix_90 = self.ipixs_in_percentage(self.skymap,0.9)
                self.ra, self.dec, self.extension = None, None, None
                self.source_type = "Prior"
            except:
                message = 'Location must be one of the following)\n'
                message += '\t(1) RA, DEC in degrees\n'
                message += '\t(2) URL to Skymap .fits file\n'
                message += '\t(3) Path to Skymap on cobalt'
                print(message)
                sys.exit()

        start = Time(dateutil.parser.parse(tstart)).mjd
        stop = Time(dateutil.parser.parse(tstop)).mjd
        
        start_date = Time(dateutil.parser.parse(tstart)).datetime
        start_str = '{:02d}'.format(start_date.year) + '_' + '{:02d}'.format(start_date.month) + '_' + '{:02d}'.format(start_date.day)
        
        self.start = start
        self.stop = stop
        self.duration = stop - start
        self.centertime = (start + stop) / 2.
        self.trigger = kwargs.pop("trigger", self.centertime)
        self.alert_event = kwargs.pop('alert_event', False)

        dirname = os.environ.get('FAST_RESPONSE_OUTPUT')
        if dirname is None:
            self.dirname = os.getcwd()
        else:
            self.dirname = dirname
        #self.dirname = kwargs.get('dirname', os.getcwd())
        self.analysisid = start_str + '_' + self.name.replace(' ', '_') 
        self.analysispath = self.dirname + '/' + self.analysisid

        self.save_output = kwargs.pop("save", True)

        if os.path.isdir(self.analysispath) and self.save_output:
            print("Directory {} already exists. Deleting it ...".format(self.analysispath))
            subprocess.call(['rm', '-r', self.analysispath])
            subprocess.call(['mkdir',self.analysispath])
            #sys.exit()
        elif self.save_output:
            subprocess.call(['mkdir',self.analysispath])

        skip_events = kwargs.pop("Skipped Events", None)
        self.skipped = skip_events
        self.skipped_event = None
        if 'test' in self.name.lower():
            self.scramble = True
            print('Working on scrambled data')
        else:
            self.scramble = False
            print('Working on unscrambled (UNBLINDED) data')

        print(self.intro_message())

        self.llh = self.initialize_llh(skipped = skip_events, extension = self.extension, 
                        scramble=self.scramble, alert_event=self.alert_event)
        self.ts, self.ns, self.p, self.sigma = None, None, None, None
        self.tsd = None
        self.skymap_fit_ra = None
        self.skymap_fit_dec = None
        self.inj = None
        self.coincident_events = None
        self.upperlimit, self.upperlimit_ninj = None, None
        self.ns_profile = None
        self.low5, self.high5 = None, None

  
    def initialize_llh(self, skipped = None, extension = None, scramble=False, alert_event=False):
        r''' Initializes a point source llh

        Parameters:
        -----------
        skipped: str
            Events to remove from analysis in format <RUN>:<EVENT>
        extension: float
            Source extension in radians
        Returns:
        --------
        llh: Skylab ps_llh object
        ''' 
        print("Initializing Point Source LLH in Skylab")
        start, stop = self.start, self.stop

        ###################### BEGIN DATASET ######################
        t0 = time.time()

        dataset = Datasets['GFUOnline']

        # only ends up in the if part for archival analyses, where it takes
        #a while to query i3live
        # CHANGE THIS FROM LIVESTREAM TO JUST GRAB THE SEASONS FROM THE DATASET CLASS
        if stop < 58850.:
            print("Old times, just grabbing archival data")
            exps, grls = [], []
            for season in ["IC86, 2016", "IC86, 2017", "IC86, 2018"]:
                exp, mc, livetime = dataset.season(season, floor=np.radians(0.2))
                grl = dataset.grl(season)
                exps.append(exp)
                grls.append(grl)
            exp = np.concatenate(exps)
            exp.sort(order='time')
            grl = np.concatenate(grls)
            grl.sort(order='run')
            livetime = grl['livetime'].sum()
        else:
            #print("querying the i3live database")
            exp, mc, livetime, grl = dataset.livestream(start - 6., stop,
                                                    append=["IC86, 2017", "IC86, 2018"], 
                                                    floor=np.radians(0.2))  #TEMPORARY
                                                    #THIS IS PS STANDARD, NOT FR STANDARD

        sinDec_bins = dataset.sinDec_bins("livestream")
        energy_bins = dataset.energy_bins("livestream")
        ###################### END DATASET   ######################

        ##################### BEGIN LIKELIHOOD MODELS #####################
        gamma_llh = 2.5 if alert_event else 2.0
        #print("working on energy term")
        ttt = time.time()
        llh_model = EnergyLLH(twodim_bins=[energy_bins, sinDec_bins],   # energy and sin(dec) binnings
                            allow_empty=True,                           # allow empty bins.
                            spectrum = PowerLaw(A=1, gamma=gamma_llh, E0=1000.))                               
        #print("energy term took {} seconds".format(time.time() - ttt))
        #print("working on temporal part")
        ttt = time.time()
        box = TemporalModel(grl=grl,
                            poisson_llh=True,   # set to True for GRBLLH style likelihood with poisson term
                            days=10,            # use 10 days on either side of ontime for background calc
                            signal=BoxProfile(start, stop))
        #print("temporal part took {} seconds".format(time.time() - ttt))
        if skipped is not None:
            try:
                event = skipped[0]
                mjd_keys = exp['time'][(exp['run'] == int(event[0])) & (exp['event'] == int(event[1]))]
                self.skipped_event = exp[exp['time'] == mjd_keys[0]][0]
                exp = dataset.remove_ev(exp, mjd_keys=mjd_keys[0])
            except:
                print("Was not able to remove event {}".format(skipped))
                exp = exp

        if extension is not None:
            src_extension = float(extension)
        else:
            src_extension = None

        #print("assembling it all")
        ttt = time.time()
        llh = PointSourceLLH(exp,                   # array with data events
                                mc,                    # array with Monte Carlo events
                                livetime,              # total livetime of the data events
                                ncpu=5,               # use 10 CPUs when computing trials
                                scramble=scramble,        # use scrambled data, set to False for unblinding
                                timescramble=True,     # use full time scrambling, not just RA scrambling
                                llh_model=llh_model,   # likelihood model
                                temporal_model=box,    # use box profile for temporal model
                                nsource_bounds=(0., 1e3),  # bounds on fitted number of signal events
                                src_extension = src_extension, # Model symmetrically extended source
                                nsource=1.)            # seed for nsignal fit

        #print("assembling took {} seconds".format(time.time() - ttt))                                               
        print("Initializing Point Source COMPLETE")
        print("LLH Initialization took {} seconds\n\n".format(time.time() - t0))
        return llh

    def initialize_injector(self, gamma=2.0, seed=123):
        r'''Method to make relevant injector in Skylab, done for analysis
        checks as well as for calculating sensitivities

        Parameters:
        -----------
        self
        gamma: float
            spectral index to inject
        Returns:
        --------
        inj: Skylab injector object'''
        if self.skymap is not None:
            from skylab.ps_injector import PriorInjector
            print("Initializing Prior Injector")
            if self.stop - self.start > 1.:
                spatial_prior = SpatialPrior(self.skymap, containment = 0.90)
            else:
                spatial_prior = SpatialPrior(self.skymap, containment = 0.90)
            self.spatial_prior = spatial_prior
            inj = PriorInjector(spatial_prior, gamma=gamma, e_range = (0,np.inf), 
                                    E0=1000., seed = seed) #1000 for 1 TeV
            inj.fill(self.llh.exp, self.llh.mc, self.llh.livetime, temporal_model=self.llh.temporal_model)
            self.inj = inj
            return inj
        else:
            print("Initializing Point Source Injector")
            inj = PointSourceInjector(gamma = gamma, E0 = 1000.) # injected spectrum dN/dE = A * (E / 1 TeV)^-2
            inj.fill(self.dec, self.llh.exp, self.llh.mc, self.llh.livetime, temporal_model=self.llh.temporal_model)
            self.inj = inj
            return inj

    def unblind_TS(self):
        r''' Unblind TS, either sky scan for spatial prior,
        or just at one location for a point source

        Parameters:
        -----------
        self 
        --------
        ts, ns: Test statistic and best fit ns
        ''' 
        if self.skymap is not None:
            t0 = time.time()
            if 'ice' in self.name.lower():
                spatial_prior = SpatialPrior(self.skymap, containment = 0.90)
            else:
                spatial_prior = SpatialPrior(self.skymap, containment = 0.90)
            pixels = np.arange(len(self.skymap))
            t1 = time.time()
            print("Starting scan")
            val = self.llh.scan(0.0,0.0, scramble = False,spatial_prior=spatial_prior,
                                time_mask = [self.duration/2., self.centertime],
                                pixel_scan=[self.nside, 3.0])
            t2 = time.time()
            print("finished scan, took {} s".format(t2-t1))
            try:
                ts = val['TS_spatial_prior_0'].max()
                max_prior = np.argmax(val['TS_spatial_prior_0'])
                ns = val['nsignal'][max_prior]
                if ts > 0:
                    self.skymap_fit_ra = val['ra'][max_prior]
                    self.skymap_fit_dec = val['dec'][max_prior]
                else:
                    max_pix = np.argmax(self.skymap)
                    theta, phi = hp.pix2ang(self.nside, max_pix)
                    dec = np.pi/2. - theta
                    self.skymap_fit_ra = phi
                    self.skymap_fit_dec = dec
            except:
                ts, ns = 0., 0.
                max_pix = np.argmax(self.skymap)
                theta, phi = hp.pix2ang(self.nside, max_pix)
                dec = np.pi/2. - theta
                self.skymap_fit_ra = phi
                self.skymap_fit_dec = dec

        else:
            ts, ns = self.llh.fit_source(src_ra=self.ra, src_dec=self.dec)
            spatial_weights = self.llh.llh_model.signal(self.ra, self.dec, 
                    self.llh._events, src_extension=self.extension)[0] / self.llh._events['B']
            params = ns.copy()
            params.pop('nsignal')
            energy_ratio, _ = self.llh.llh_model.weight(self.llh._events, **params)
            #print(spatial_weights)
            #print(energy_ratio)
            #print(np.max(spatial_weights * energy_ratio))
            ns = ns['nsignal']
            msk = spatial_weights * energy_ratio > 10
            if len(spatial_weights[msk]) > 0:
                self.coincident_events = []
                #print self.llh._events.dtype.names
                for ev, s_w, en_w in zip(self.llh._events[msk], 
                            spatial_weights[msk], energy_ratio[msk]):
                    self.coincident_events.append({})
                    for key in ['run', 'event', 'ra', 'dec', 'sigma', 'logE', 'time']:
                        self.coincident_events[-1][key] = ev[key]
                    del_psi = deltaPsi([ev['dec']], [ev['ra']], [self.dec], [self.ra])[0]
                    self.coincident_events[-1]['delta_psi'] = del_psi 
                    self.coincident_events[-1]['spatial_w'] = s_w
                    self.coincident_events[-1]['energy_w'] = en_w
            #print self.coincident_events
        print("TS = {}".format(ts))
        print("ns = {}\n\n".format(ns))
        self.ts, self.ns = ts, ns
        return ts, ns
   
    def calc_pvalue(self, ntrials=1000):
        r''' Given an unblinded TS value, calculate the p value

        Parameters:
        -----------
        None
        Returns:
        --------
        p: float
            P value for this analysis
        ''' 

        if self.ts is None:
            print("Need TS value to find p value")
            p = np.nan
            self.tsd = None
        elif self.ts == 0.0:
            print("TS=0, no need to run background trials")
            p = 1.0
            self.tsd = None
        else:
            print("Need to reinitialize LLH for background trials")
            self.llh = self.initialize_llh(extension = self.extension, scramble=True)
            if self.source_type == "PS":
                print("Running {} Background trials".format(ntrials))
                t0 = time.time()
                trials = self.llh.do_trials(int(ntrials), src_ra=self.ra, src_dec=self.dec)
                p = sum(trials['TS'] >= self.ts) / float(trials.size)
                self.tsd = trials['TS']
                t1 = time.time()
                print("Background Trials complete, p = {}".format(p))
                print("Background trials took {} seconds ({} s / trial)".format(t1 - t0, (t1-t0) / float(ntrials)))
            else:
                print("Running {} Background trials".format(ntrials))
                t0 = time.time()
                tsd = self.run_prior_bg(ntrials)
                self.tsd = tsd
                p = sum(tsd >= self.ts) / float(tsd.size)
                t1 = time.time()
                print("Background Trials complete, p = {}".format(p))
                print("Background trials took {} seconds ({} s / trial)".format(t1 - t0, (t1-t0) / float(ntrials)))

        self.p = p
        sigma = self.significance()
        print("Significance: {:.3f}sigma\n\n".format(sigma))
        return p
                
    def run_prior_bg(self, ntrials = 1000):
        r''' helper method to run background trials for the spatial prior analyses

        Returns:
        --------
        tsd: array-like
            TS distribution for the background only trials
        ''' 
        tsd = []
        spatial_prior = SpatialPrior(self.skymap, containment=0.90)
        toolbar_width = 50
        sys.stdout.write("[%s]" % (" " * toolbar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (toolbar_width+1))
        t = time.time()
        for i in range(int(ntrials)):
            if ntrials > 50 and i % (ntrials / 50) == 0:
                sys.stdout.write("#")
                sys.stdout.flush()
            val = self.llh.scan(0.0, 0.0, scramble = True, seed=i, spatial_prior=spatial_prior,
                            time_mask = [self.duration/2., self.centertime], pixel_scan=[self.nside,3.0])
            try:
                tsd.append(val['TS_spatial_prior_0'].max())
            except ValueError:
                tsd.append(0.)
        tsd = np.asarray(tsd,dtype=float)
        t1 = time.time()
        print("{} trials took {}s".format(ntrials, t1-t))
        return tsd

    def significance(self):
        r'''Given p value, report significance

        Parameters:
        ----------
        p: float
            P value
        Returns:
        --------
        sigma: float
            Significance
        '''
        if self.p is not None:
            sigma = np.sqrt(2)*erfinv(1-2*self.p)
            self.sigma = sigma
            return sigma
        else:
            print("Need a p value to calculate significance")
            return np.nan

    def ns_scan(self, params = {'spectrum': 'dN/dE = 1.00e+00 * (E / 1.00e+03 GeV)^-2.00 [GeV^-1cm^-2s^-1]'}):
        r''' Calculate the llh as a function of ns, to 
        see how shallow the likelihood space is.

        Parameters:
        -----------
        params: mappable
            params dictionary for skylab.ps_llh object. Default is 
            that of this analysis
        Returns:
        ---------
        xs, llhs: array-like
            Values of ns scanned and corresponding llh values'''
        bounds = self.ns * 2.5
        if bounds == 0.0:
            bounds = 1.0
        xs = np.linspace(0., bounds, 101)
        llhs = []
        for n in xs:
            llhs.append(self.llh.llh(n, **params)[0])
        llhs = np.array(llhs)

        if self.save_output:
            fig, ax = plt.subplots()
            plt.plot(xs, -llhs, lw = 3)
            plt.axvline(self.ns, ls = '--', color = 'gray', lw = 1)
            plt.xlabel(r'$n_s$', fontsize = 18)
            plt.ylabel(r'$-\log \frac{\mathcal{L}}{\mathcal{L}_0}$', fontsize = 18)
            plt.savefig(self.analysispath + '/llh_ns_scan.png', bbox_inches='tight', dpi=200)

        self.ns_profile = (xs, llhs)
        return xs, llhs

    def upper_limit(self, n_per_sig = 500, p0 = None):
        r'''After calculating TS, find upper limit
        Assuming an E^-2 spectrum
        Returns:
        --------
        Value of E^2 dN / dE in units of TeV / cm^2 
        '''
        print("Beginning upper limit calculation")
        if self.inj is None:
            gamma_inj = 2.5 if self.alert_event else 2.0
            self.initialize_injector(gamma = gamma_inj)
        if self.skymap is None:
            #ninj = np.linspace(1., 6., 6)
            ninj = np.array([1., 1.5, 2., 2.5, 3., 4., 5., 6.])
            passing = []
            for n in ninj:
                results = self.llh.do_trials(n_per_sig, src_ra=self.ra, src_dec=self.dec,
                            injector=self.inj, mean_signal=n, poisson=True)
                msk = results['TS'] > self.ts
                npass = len(results['TS'][msk])
                passing.append((n, npass, n_per_sig))
        else:
            print('Upper limit with spatial prior not yet implemented')
            return
        signal_fluxes, passing, number = zip(*passing)
        signal_fluxe = np.array(signal_fluxes)
        passing = np.array(passing, dtype=float)
        number = np.array(number, dtype=float)
        passing /= number
        errs = utils.binomial_error(passing, number)
        fits, plist = [], []
        try:
            fits.append(sensitivity_fit(signal_fluxes, passing, errs, utils.chi2cdf, p0=p0))
            plist.append(fits[-1]['pval'])
        except:
            print("chi2 fit failed in upper limit calculation")
        try:
            fits.append(sensitivity_fit(signal_fluxes, passing, errs, utils.erfunc, p0=p0))
            plist.append(fits[-1]['pval'])
        except:
            print("Error function fit failed in upper limit calculation")
        try:
            fits.append(sensitivity_fit(signal_fluxes, passing, errs, utils.incomplete_gamma, p0=p0))
            plist.append(fits[-1]['pval'])
        except:
            print("Incomplte gamma fit failed in upper limit calculation")
        #Find best fit of the three, make it look different in plot
        plist = np.array(plist)
        best_fit_ind= np.argmax(plist)
        fits[best_fit_ind]['ls'] = '-'
        self.upperlimit = self.inj.mu2flux(fits[best_fit_ind]['sens'])
        self.upperlimit_ninj = fits[best_fit_ind]['sens']

        fig, ax = plt.subplots()
        for fit_dict in fits:
            ax.plot(fit_dict['xfit'], fit_dict['yfit'], 
                     label = r'{}: $\chi^2$ = {:.2f}, d.o.f. = {}'.format(fit_dict['name'], fit_dict['chi2'], fit_dict['dof']),
                    ls = fit_dict['ls'])
            if fit_dict['ls'] == '-':
                ax.axhline(0.9, color = 'm', linewidth = 0.3, linestyle = '-.')
                ax.axvline(fit_dict['sens'], color = 'm', linewidth = 0.3, linestyle = '-.')
                ax.text(3.5, 0.8, r'Sens. = {:.2f} events'.format(fit_dict['sens']), fontsize = 16)
                ax.text(3.5, 0.7, r' = {:.1e}'.format(self.upperlimit * self.duration * 86400. * 1e6) + r' GeV cm$^{-2}$', fontsize = 16)
                #ax.text(5, 0.5, r'Sens. = {:.2e}'.format(self.inj.mu2flux(fit_dict['sens'])) + ' GeV^-1 cm^-2 s^-1')
        ax.errorbar(signal_fluxes, passing, yerr=errs, capsize = 3, linestyle='', marker = 's', markersize = 2)
        ax.legend(loc=4, fontsize = 14)
        ax.set_xlabel(r'$\langle n_{inj} \rangle$', fontsize = 14)
        ax.set_ylabel(r'Fraction TS $>$ unblinded TS', fontsize = 14)
        if self.save_output:
            plt.savefig(self.analysispath + '/upper_limit_distribution.png', bbox_inches='tight', dpi=200)
        print("Found upper limit of {:.2f} events".format(self.upperlimit_ninj))
        return self.upperlimit

    """def write_circular(self):
        r'''After running analysis and calculating upper
            limit, write a GCN circular'''
        ############################FINISH THIS PART ######################
        with open('./circular_shell.txt', 'r') as f:
            circular_text = f.readlines()
            f.close() 
        new_circ = []
        for line in circular_text:
            if '' in line:
                new_l = line.replace('', '')
        
        with open(self.analysis_path + '/' + self.analysisid + 'gcn_circular.txt', 'w') as f:
            for line in new_circ:
                f.write(line)"""

    def intro_message(self):
        r'''Print a message with info about the source once
        the analysis is running'''
        int_str = '*'*80
        int_str += '\n*' + ' '*78 + '*\n'
        int_str += '*' + ' '*((78-len(self.name))/2) + self.name + ' '*((78-len(self.name))/2 + len(self.name)%2) + '*'
        int_str += '\n*' + ' '*78 + '*\n'
        int_str += '*'*80 + '\n'
        int_str += '  '*5 + str(Time(self.start, format = 'mjd', scale = 'utc').iso)
        int_str += '  '*5 + str(Time(self.stop, format = 'mjd', scale = 'utc').iso) + '\n'
        if self.skymap is not None:
            int_str += ' '*10 + 'Skymap file:' + self.skymap_url
        else:
            int_str += ' '*10 + 'RA, dec.: {:.2f}, {:+.3f}'.format(self.ra*180. / np.pi, self.dec*180. / np.pi)
            exts = 0. if self.extension is None else self.extension
            int_str += ' '*10 + 'Extension: {}'.format(exts * 180. / np.pi)
        int_str += '\n\n'
        return int_str

    def get_dirname(self):
        r'''Returns analysis directory name
        Returns:
        --------
        dirname: str
            Analysis directory name
        '''
        return self.analysisid

    def save_results(self):
        r'''Once analysis has been performed, push
        all relevant info to a new directory
        Parameters:
        -----------
        None
        Returns:
        ________
        None
        '''
        results = {}
        for key, val in [('start', self.start), ('stop', self.stop), 
                        ('name', self.name),  ('source_type', self.source_type),
                        ('analysisid', self.analysisid), ('analysispath', self.analysispath)]:
            results[key] = val
        for key, val in [('ts', self.ts), ('p', self.p), ('sigma', self.sigma), ('ns', self.ns),
                        ('tsd', self.tsd), ('extension', self.extension), 
                        ('skymap', self.skymap), ('ra', self.ra), ('dec', self.dec),
                        ('coincident_events', self.coincident_events), ('skipped', self.skipped_event), 
                        ('upper_limit', self.upperlimit), ('upper_limit_ninj', self.upperlimit_ninj),
                        ('ns_profile', self.ns_profile), ('low_en', self.low5), ('high_en', self.high5)]:
            if val is not None:
                results[key] = val
        print("Saving results to directory:\n\t{}".format(self.analysispath))
        with open(self.analysispath + '/' + self.analysisid + '_' + 'results.pickle', 'wb') as f:
            pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)
        print("Results successfully saved")
        return results       
 
    def generate_report(self):
        r'''Generates report using class attributes
        and the ReportGenerator Class'''
        #report = ReportGenerator(self.name, self.trigger, self.start, self.stop,
        #                            self.ts, self.ns, self.source_type, self.analysisid, **vars(self))

        report_kwargs = vars(self).copy()
        print("\nGenerating PDF Report")
        for key in ['name', 'trigger', 'start', 'stop', 'ts', 'ns', 'source_type', 'analysisid']:
            report_kwargs.pop(key)
        report = ReportGenerator(self.name, self.trigger, self.start, self.stop, 
                                self.ts, self.ns, self.source_type, self.analysisid, **report_kwargs)
        report.generate_report()
        report.make_pdf()
        print("Report saved to output directory")
        username = subprocess.check_output(['whoami']).strip('\n')
        if os.path.isdir('/home/{}/public_html/FastResponse/'.format(username)):
            print("Copying report to {}'s public html in FastResponse Directory".format(username))
            shutil.copy2(self.analysispath + '/' + self.analysisid + "_report.pdf", 
                '/home/{}/public_html/FastResponse/{}_report.pdf'.format(username, self.analysisid))
        else:
            print("Creating FastResponse Directory in {}'s public html and copyting report")
            os.mkdir('/home/{}/public_html/FastResponse/'.format(username))
            shutil.copy2(self.analysispath + '/' + self.analysisid + "_report.pdf", 
                '/home/{}/public_html/FastResponse/{}_report.pdf'.format(username, self.analysisid))

    def plot_tsd(self):
        r'''Outputs a plot of the background TS distributions
        as well as the observed TS
        '''
        if self.tsd is not None:
            fig, ax = plt.subplots()
            plt.hist(self.tsd, bins= np.linspace(0., 25, 30), 
                    label="Background Scrambles", density=True)
            plt.axvline(self.ts, color = 'k', label = "Observed TS")
            plt.legend(loc=1)
            plt.xlabel("TS")
            plt.ylabel("Probability Density")
            plt.yscale('log')
            plt.savefig(self.analysispath + '/TS_distribution.png', bbox_inches='tight')

    def plot_ontime(self):
        r'''Plots ontime events with either the full sky spatial
        prior or a zoomed in version like traditional fast response
        followups
        '''
        self.plot_skymap_zoom()
        self.plot_skymap()

    def plot_skymap_zoom(self):
        r'''Make a zoomed in portion of a skymap with
        all neutrino events within a certain range
        '''
        events = self.llh.exp
        events = events[(events['time'] < self.stop) & (events['time'] > self.start)]
        #print(np.count_nonzero(events['sigma'] * 180. / np.pi < 0.2))

        col_num = 5000
        seq_palette = sns.color_palette("Spectral", col_num)
        lscmap = mpl.colors.ListedColormap(seq_palette)

        rel_t = np.array((events['time'] - self.start) * col_num / (self.stop - self.start), dtype = int)
        cols = np.array([seq_palette[j] for j in rel_t])

        if self.skymap is not None:
            skymap = self.skymap
            ra = self.skymap_fit_ra
            dec = self.skymap_fit_dec
            if 'ice' in self.name.lower():
                skymap = np.log10(skymap)
            label_str = "Scan Hot Spot"
            cmap = None
        else:
            skymap = np.zeros(hp.nside2npix(256))
            ra = self.ra
            dec = self.dec
            label_str = self.name
            cmap = mpl.colors.ListedColormap([(1.,1.,1.)] * 50)

        plot_zoom(skymap, ra, dec, "Zoomed Scan Results", range = [0,10], reso=3., cmap = cmap)

        if self.skipped is not None:
            try:
                msk = events['run'] == int(self.skipped[0][0])
                msk *= events['event'] == int(self.skipped[0][1])
                plot_events(self.skipped_event['dec'], self.skipped_event['ra'], self.skipped_event['sigma']*2.145966, 
                            ra, dec, 2*6, sigma_scale=1.0, constant_sigma=False, same_marker=True, 
                            energy_size=True, col = 'k', with_dash=True)
                events = events[~msk]
                cols = cols[~msk]
            except:
                print("Removed event not in GFU")

        if (self.stop - self.start) <= 14.:
            plot_events(events['dec'], events['ra'], events['sigma']*2.145966, ra, dec, 2*6, sigma_scale=1.0,
                    constant_sigma=False, same_marker=True, energy_size=True, col = cols)
                    #2.145966 for 90% containment
        else:
            #Long time windows means don't plot contours
            plot_events(events['dec'], events['ra'], events['sigma']*2.145966, ra, dec, 2*6, sigma_scale=None,
                    constant_sigma=False, same_marker=True, energy_size=True, col = cols)

        plt.scatter(0,0, marker='*', c = 'k', s = 130, label = label_str) 
        #hp.projscatter(np.pi/2-dec, ra, marker='o', linewidth=2, edgecolor='k', linestyle=':', facecolor="None", s=5200*1, alpha=1.0)
        #hp.projscatter(np.pi/2-dec, ra, marker='o', linewidth=2, edgecolor='g', linestyle=':', facecolor="None", s=5200*4.0, alpha=1.0)

        plt.text(1.2*np.pi / 180., 2.8*np.pi / 180., 'IceCube\nPreliminary', color = 'r', fontsize = 22)
        plt.legend(loc = 2, ncol=2, mode = 'expand', fontsize = 15.5, framealpha = 0.95)
        plot_color_bar(range=[0,6], cmap=lscmap, col_label=r"IceCube Event Time",
                    offset=-55, labels = [r'-$\Delta t \Bigg/ 2$', r'+$\Delta t \Bigg/ 2$'])
        plt.savefig(self.analysispath + '/' + self.analysisid + 'unblinded_skymap_zoom.png',bbox_inches='tight')
        plt.savefig(self.analysispath + '/' + self.analysisid + 'unblinded_skymap_zoom.pdf',bbox_inches='tight', dpi=300)

    def plot_skymap(self):
        r''' Make skymap with event localization and all
        neutrino events on the sky within the given time window
        '''

        events = self.llh.exp
        events = events[(events['time'] < self.stop) & (events['time'] > self.start)]

        col_num = 5000
        seq_palette = sns.color_palette("Spectral", col_num)
        lscmap = mpl.colors.ListedColormap(seq_palette)

        rel_t = np.array((events['time'] - self.start) * col_num / (self.stop - self.start), dtype = int)
        cols = [seq_palette[j] for j in rel_t]

        # Set color map and plot skymap
        pdf_palette = sns.color_palette("Blues", 500)
        cmap = mpl.colors.ListedColormap(pdf_palette)
        cmap.set_under("w")

        if self.skymap is None:
            skymap = np.zeros(hp.nside2npix(256))
            max_val = 1.
            cmap = mpl.colors.ListedColormap([(1.,1.,1.)] * 50)
        else:
            if 'ice' in self.name.lower():
                skymap = np.log10(self.skymap)
                max_val = max(skymap)
            else:
                skymap = self.skymap
                max_val = max(skymap)
        hp.mollview(skymap,coord='C',cmap=cmap)
        hp.graticule(verbose=False)

        theta=np.pi/2 - events['dec']
        phi = events['ra']

        #plot 90% containment errors
        sigma_90 = events['sigma']*2.145966

        # plot events on sky with error contours
        hp.projscatter(theta,phi,c=cols,marker='x',label='GFU Event',coord='C')

        if (self.stop - self.start) <= 0.5:        #Only plot contours if less than 2 days
            for i in range(events['ra'].size):
                my_contour = self.contour(events['ra'][i], 
                                    events['dec'][i],sigma_90[i], 256)
                hp.projplot(my_contour[0], my_contour[1], linewidth=2., 
                                    color=cols[i], linestyle="solid",coord='C')

        if self.skymap is None:
            src_theta = np.pi/2. - self.dec
            src_phi = self.ra
            hp.projscatter(src_theta, src_phi, c = 'k', marker = '*',
                                label = self.name, coord='C', s=350)

        plt.title('Fast Response Skymap')
        plt.legend(loc=1)
        plt.savefig(self.analysispath + '/' + self.analysisid + 'unblinded_skymap.png',bbox_inches='tight')

    def contour(self,ra, dec, sigma,nside):
        r''' Function for plotting contours on skymaps

        Parameters:
        -----------
        ra: ndarray
            Array of ra for events 
        dec: ndarray
            Array of dec for events
        sigma: ndarray
            Array of sigma to make contours around events
        nside:
            nside of healpy map
        Returns:
        --------
        Theta: array
            array of theta values of contour
        Phi: array
            array of phi values of contour 
        '''
        dec = np.pi/2 - dec
        sigma = np.rad2deg(sigma)
        delta, step, bins = 0, 0, 0
        delta= sigma/180.0*np.pi
        step = 1./np.sin(delta)/20.
        bins = int(360./step)
        Theta = np.zeros(bins+1, dtype=np.double)
        Phi = np.zeros(bins+1, dtype=np.double)
        # define the contour
        for j in range(0,bins) :
                phi = j*step/180.*np.pi
                vx = np.cos(phi)*np.sin(ra)*np.sin(delta) + np.cos(ra)*(np.cos(delta)*np.sin(dec) + np.cos(dec)*np.sin(delta)*np.sin(phi))
                vy = np.cos(delta)*np.sin(dec)*np.sin(ra) + np.sin(delta)*(-np.cos(ra)*np.cos(phi) + np.cos(dec)*np.sin(ra)*np.sin(phi))
                vz = np.cos(dec)*np.cos(delta) - np.sin(dec)*np.sin(delta)*np.sin(phi)
                idx = hp.vec2pix(nside, vx, vy, vz)
                DEC, RA = hp.pix2ang(nside, idx)
                Theta[j] = DEC
                Phi[j] = RA
        Theta[bins] = Theta[0]
        Phi[bins] = Phi[0]

        return Theta, Phi

    def make_dNdE(self):
        r'''Make an E^-2 dNdE with the central 90% 
        for the most relevant declination'''
        if self.skymap is None:
            dec_mask_1 = self.llh.mc['dec'] > self.dec - (5. * np.pi / 180.)
            dec_mask_2 = self.llh.mc['dec'] < self.dec + (5. * np.pi / 180.)
            dec_mask_3, dec_mask_4 = None, None
        else:
            min_dec, max_dec = self.dec_skymap_range()
            dec_mask_1 = self.llh.mc['dec'] > min_dec - (5. * np.pi / 180.)
            dec_mask_2 = self.llh.mc['dec'] < min_dec + (5. * np.pi / 180.)
            dec_mask_3 = self.llh.mc['dec'] > max_dec - (5. * np.pi / 180.)
            dec_mask_4 = self.llh.mc['dec'] < max_dec + (5. * np.pi / 180.)
        dec_mask = dec_mask_1 * dec_mask_2
        fig, ax = plt.subplots(figsize = (8,5))
        fig.set_facecolor('white')
        lab = 'Min dec.' if self.skymap is not None else None
        a = plt.hist(self.llh.mc['trueE'][dec_mask], bins = np.logspace(1., 8., 50), 
                weights = self.llh.mc['ow'][dec_mask] * np.power(self.llh.mc['trueE'][dec_mask], -1.) / self.llh.mc['trueE'][dec_mask], 
                histtype = 'step', linewidth = 2., color = sns.xkcd_rgb['windows blue'], label = lab)
        plt.yscale('log')
        plt.xscale('log')
        plt.grid(which = 'major', alpha = 0.25)
        plt.xlabel('Energy (GeV)', fontsize = 24)
        cdf = np.cumsum(a[0]) / np.sum(a[0])
        low_5 = np.interp(0.05, cdf, a[1][:-1])
        median = np.interp(0.5, cdf, a[1][:-1])
        high_5 = np.interp(0.95, cdf, a[1][:-1])
        self.low5 = low_5
        self.high5 = high_5
        plt.axvspan(low_5, high_5, color = sns.xkcd_rgb['windows blue'], alpha = 0.25, label="Central 90\%")
        lab = 'Median (min dec.)' if self.skymap is not None else 'Median'
        plt.axvline(median, c = sns.xkcd_rgb['windows blue'], alpha = 0.75, label = lab)
        if dec_mask_3 is not None:
            dec_mask = dec_mask_3 * dec_mask_4
            a = plt.hist(self.llh.mc['trueE'][dec_mask], bins = np.logspace(1., 8., 50), 
                weights = self.llh.mc['ow'][dec_mask] * np.power(self.llh.mc['trueE'][dec_mask], -1.) / self.llh.mc['trueE'][dec_mask], 
                histtype = 'step', linewidth = 2., color = sns.xkcd_rgb['dark navy blue'], label = 'Max dec.')
            cdf = np.cumsum(a[0]) / np.sum(a[0])
            low_5 = np.interp(0.05, cdf, a[1][:-1])
            median = np.interp(0.5, cdf, a[1][:-1])
            high_5 = np.interp(0.95, cdf, a[1][:-1])
            plt.axvspan(low_5, high_5, color = sns.xkcd_rgb['dark navy blue'], alpha = 0.25)
            plt.axvline(median, c = sns.xkcd_rgb['dark navy blue'], alpha = 0.75, label = "Median (max dec.)", ls = '--')
        plt.xlim(1e1, 1e8)
        plt.legend(loc=4, fontsize=18)
        plt.savefig(self.analysispath + '/central_90_dNdE.png',bbox_inches='tight')
        

    def ipixs_in_percentage(self,skymap,percentage):
        """Finding ipix indices confined in a given percentage.
        
        Input parameters
        ----------------
        skymap: ndarray
            array of probabilities for a skymap
        percentage : float
            fractional percentage from 0 to 1  
        Return
        ------- 
        ipix : numpy array
            indices of pixels within percentage containment
        """
    
        sort = sorted(skymap, reverse = True)
        cumsum = np.cumsum(sort)
        index, value = min(enumerate(cumsum),key=lambda x:abs(x[1]-percentage))

        index_hpx = range(0,len(skymap))
        hpx_index = np.c_[skymap, index_hpx]

        sort_2array = sorted(hpx_index,key=lambda x: x[0],reverse=True)
        value_contour = sort_2array[0:index]

        j = 1 
        table_ipix_contour = [ ]
        for i in range (0, len(value_contour)):
            ipix_contour = int(value_contour[i][j])
            table_ipix_contour.append(ipix_contour)
    
        ipix = table_ipix_contour
          
        return np.asarray(ipix,dtype=int)

    def dec_skymap_range(self):
        r''' Compute minimum and maximum declinations within
        of the 90% contour of a given skymap

        Returns:
        --------
        low: float
            lowest declination in 90%
        high: float
            highest declination in 90%
        '''

        src_theta, src_phi = hp.pix2ang(self.nside,self.ipix_90)
        src_dec = np.pi/2. - src_theta
        src_dec = np.unique(src_dec)

        low = src_dec.min()
        high = src_dec.max()

        return low, high


def deltaPsi(dec1, ra1, dec2, ra2):
    """
    Calculate angular distance.
    
    Args:
        dec1: Declination of first direction in radian
        ra1: Right ascension of first direction in radian
        dec2: Declination of second direction in radian
        ra2: Right ascension of second direction in radian
        
    Returns angular distance in radian
    """
    return deltaPsi2(np.sin(dec1), np.cos(dec1), np.sin(ra1), np.cos(ra1), np.sin(dec2), np.cos(dec2), np.sin(ra2), np.cos(ra2))

def deltaPsi2(sDec1, cDec1, sRa1, cRa1, sDec2, cDec2, sRa2, cRa2):
    """
    Calculate angular distance.
    
    Args:
        sDec1: sin(Declination of first direction)
        cDec1: cos(Declination of first direction)
        sRa1: sin(Right ascension of first direction)
        cRa1: cos(Right ascension of first direction)
        sDec2: sin(Declination of second direction)
        cDec2: cos(Declination of second direction)
        sRa2: sin(Right ascension of second direction)
        cRa2: cos(Right ascension of second direction)
        
    Returns angular distance in radian
    """
    tmp = cDec1*cRa1*cDec2*cRa2 + cDec1*sRa1*cDec2*sRa2 + sDec1*sDec2
    tmp[tmp>1.] = 1.
    tmp[tmp<-1.] = -1.
    return np.arccos(tmp)
                                 
def plot_zoom(scan, ra, dec, title, reso=3, var="pVal", range=[0, 6],cmap=None):
    if cmap is None:
        pdf_palette = sns.color_palette("Blues", 500)
        cmap = mpl.colors.ListedColormap(pdf_palette)
    hp.gnomview(scan, rot=(np.degrees(ra), np.degrees(dec), 0),
                    cmap=cmap,
                    max=max(scan)*0.1,
                    reso=reso,
                    title=title,
                    notext=True,
                    cbar=False
                    #unit=r""
                    )

    plt.plot(4.95/3.*reso*np.radians([-1, 1, 1, -1, -1]), 4.95/3.*reso*np.radians([1, 1, -1, -1, 1]), color="k", ls="-", lw=3)
    hp.graticule(verbose=False)
    plot_labels(dec, ra, reso)

def plot_color_bar(labels=[0.,2.,4.,6.], col_label=r"IceCube Event Time", range=[0,6], cmap=None, offset=-35):
    fig = plt.gcf()
    #ax = fig.add_axes([0.25, -0.03, 0.5, 0.03])
    ax = fig.add_axes([0.95, 0.2, 0.03, 0.6])
    labels = labels
    cb = mpl.colorbar.ColorbarBase(ax, cmap=ps_map if cmap is None else cmap,
                        #norm=mpl.colors.Normalize(vmin=range[0], vmax=range[1]), 
                        orientation="vertical")
    #cb.ax.minorticks_on()

    cb.set_label(col_label, labelpad=offset)
    cb.set_ticks([0., 1.])
    cb.set_ticklabels(labels)
    cb.update_ticks()
    #cb.ax.get_xaxis().set_ticklabels(labels)

def plot_labels(src_dec, src_ra, reso):
    """Add labels to healpy zoom"""
    fontsize = 20
    plt.text(-1*np.radians(1.75*reso),np.radians(0), r"%.2f$^{\circ}$"%(np.degrees(src_dec)),
             horizontalalignment='right',
             verticalalignment='center', fontsize=fontsize)
    plt.text(-1*np.radians(1.75*reso),np.radians(reso), r"%.2f$^{\circ}$"%(reso+np.degrees(src_dec)),
             horizontalalignment='right',
             verticalalignment='center', fontsize=fontsize)
    plt.text(-1*np.radians(1.75*reso),np.radians(-reso), r"%.2f$^{\circ}$"%(-reso+np.degrees(src_dec)),
             horizontalalignment='right',
             verticalalignment='center', fontsize=fontsize)
    plt.text(np.radians(0),np.radians(-1.75*reso), r"%.2f$^{\circ}$"%(np.degrees(src_ra)),
             horizontalalignment='center',
             verticalalignment='top', fontsize=fontsize)
    plt.text(np.radians(reso),np.radians(-1.75*reso), r"%.2f$^{\circ}$"%(-reso+np.degrees(src_ra)),
             horizontalalignment='center',
             verticalalignment='top', fontsize=fontsize)
    plt.text(np.radians(-reso),np.radians(-1.75*reso), r"%.2f$^{\circ}$"%(reso+np.degrees(src_ra)),
             horizontalalignment='center',
             verticalalignment='top', fontsize=fontsize)
    plt.text(-1*np.radians(2.3*reso), np.radians(0), r"declination", 
                ha='center', va='center', rotation=90, fontsize=fontsize)
    plt.text(np.radians(0), np.radians(-2*reso), r"right ascension", 
                ha='center', va='center', fontsize=fontsize)

def plot_events(dec, ra, sigmas, src_ra, src_dec, reso, sigma_scale=5., col = 'k', constant_sigma=False,
                    same_marker=False, energy_size=False, with_mark=True, with_dash=False):
    """Adds events to a healpy zoom, get events from llh."""
    cos_ev = np.cos(dec)
    tmp = np.cos(src_ra - ra) * np.cos(src_dec) * cos_ev + np.sin(src_dec) * np.sin(dec)
    dist = np.arccos(tmp)

    if sigma_scale is not None:
        sigma = np.degrees(sigmas)/sigma_scale
        sizes = 5200*sigma**2
        if constant_sigma:
            sizes = 20*np.ones_like(sizes)
        if with_dash:
            hp.projscatter(np.pi/2-dec, ra, marker='o', linewidth=2, edgecolor=col, linestyle=':', facecolor="None", s=sizes, alpha=1.0)
        else:
            hp.projscatter(np.pi/2-dec, ra, marker='o', linewidth=2, edgecolor=col, facecolor="None", s=sizes, alpha=1.0)
    if with_mark:
        hp.projscatter(np.pi/2-dec, ra, marker='x', linewidth=2, edgecolor=col, facecolor=col, s=60, alpha=1.0)

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def sensitivity_curve(ns, passing, in_ns = True, month=6, with_err = True, trim=-1, 
                            ax = None, p0 = None, fontsize = 16):
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
    
def calc_sensitivity(deltaT, sinDec, gamma, threshold = 0.5, in_ns = True, month=6, with_err = True, trim=-1, p0=None):
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
    #Find best fit of the three, make it look different in plot
    plist = np.array(plist)
    best_fit_ind= np.argmax(plist)
    return fits[best_fit_ind]
    
def sensitivity_fit(signal_fluxes, passing, errs, fit_func, p0 = None, conf_lev = 0.9):
    r'''
    Given an array of injected fluxes (or events) and passing fraction
    relative to the unblinded TS, do a fit given a CDF
    '''
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
    sens = xfit[find_nearest_idx(yfit, 0.9)]
    return {'popt': popt, 'pcov': pcov, 'chi2': chi2, 
            'dof': dof, 'xfit': xfit, 'yfit': yfit, 
            'name': name, 'pval':pval, 'ls':'--', 'sens': sens}

def pvals_for_signal(deltaT, sinDec, gamma, ns = 1, month=6, sigma = False):
    bg_trials = np.load('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/analysis_checks/bg/ps_sinDec_{}_deltaT_{}_month_{}.npy'.format(sinDec, deltaT, month))
    signal_trials = np.load('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/analysis_checks/sensitivity/nsignal_sinDec_{}_deltaT_{}_month_{}.npy'.format(sinDec, deltaT, month))
    signal_trials = signal_trials[signal_trials['gamma'] == gamma]
    signal_trials = signal_trials[signal_trials['n_inj'] == ns]
    #print(len(bg_trials['TS']))
    pvals = [100. - sp.stats.percentileofscore(bg_trials, ts, kind='strict') for ts in signal_trials['TS']]
    pvals = np.array(pvals)*0.01
    pvals = np.where(pvals==0, 1e-6, pvals)
    if not sigma:
        return pvals
    else:
        return sp.stats.norm.ppf(1. - (pvals / 2.))

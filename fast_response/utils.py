import numpy as np
import scipy as sp
from scipy.optimize       import curve_fit
from scipy.stats          import chi2
import pandas as pd
import subprocess
import os
from astropy.time import Time
import datetime
import matplotlib as mpl
mpl.use('agg')
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

############################# Plotting Parameters #############################
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
############################# Plotting Parameters #############################


def updateFastResponseWeb(analysis):
    r'''
    Create analysis specific page, and update
    plots with information from all analyses
    Parameters:
    -----------
    analysis: FastResponseAnalysis instance
    '''
    updateDataFrame(analysis)
    createFastResponsePage(analysis)
    updateFastResponseTable(analysis)
    updateFastResponsePlots()
    sync_to_roc()

def updateDataFrame(analysis):
    r'''
    Read in official Fast Response Dataframe,
    add these results, save'''
    df = pd.read_pickle('/data/user/apizzuto/fast_response_skylab/results_dataframe.pkl')
    evid = None if 'skipped' not in analysis.keys() else str(analysis['skipped']['run']) + ':' + str(analysis['skipped']['event'])
    dec = np.nan if 'dec' not in analysis.keys() else analysis['dec'] * 180. / np.pi
    ra = np.nan if 'ra' not in analysis.keys() else analysis['ra'] * 180. / np.pi
    ext = 0.0 if 'extension' not in analysis.keys() else analysis['extension'] * 180. / np.pi
    upper_lim = np.nan if 'upper_limit' not in analysis.keys() else analysis['upper_limit']
    low_en = np.nan if 'low_en' not in analysis.keys() else analysis['low_en']
    high_en = np.nan if 'high_en' not in analysis.keys() else analysis['high_en']
    new_list = [ra, dec, analysis['p'],
            analysis['ns'], pd.Timestamp(Time(analysis['start'], format='mjd').iso), 
            pd.Timedelta(analysis['stop'] - analysis['start'], unit='day'),
            ext, None, analysis['ts'], evid, upper_lim, (low_en, high_en)]
    if analysis['name'] in df.index:
        num = np.count_nonzero(df.index == analysis['name'])
        analysis['name'] += '_{}'.format(num)
    df.loc[analysis['name']] = new_list
    df.to_pickle('/data/user/apizzuto/fast_response_skylab/results_dataframe.pkl')    

def createFastResponsePage(analysis):
    r'''
    Create analysis specific page
    Parameters:
    -----------
    analysis: FastResponseAnalysis results from pickle file
    '''
    new_f = []
    keypairs = [('ANALYSISTS', 'ts'), ('ANALYSISNS', 'ns'), ('ANALYSISP', 'p')]
    with open('/data/user/apizzuto/fast_response_skylab/fast-response/fast_response/html/analysis_base.html', 'r') as f:
        for line in f.readlines():
            for k, r in keypairs:
                if k in line:
                    line = line.replace(k, '{:.3f}'.format(analysis[r]))
            new_f.append(line)
        for i in range(len(new_f)):
            if 'ANALYSISCREATIONTIME' in new_f[i]:
                new_f[i] = new_f[i].replace('ANALYSISCREATIONTIME', str(datetime.datetime.utcnow())[:-7])
            if 'ANALYSISSTART' in new_f[i]:
                new_f[i] = new_f[i].replace('ANALYSISSTART', Time(analysis['start'], format='mjd').iso)
            if 'ANALYSISDURATION' in new_f[i]:
                dur = (analysis['stop'] - analysis['start']) * 86400.
                new_f[i] = new_f[i].replace('ANALYSISDURATION', str(dur))
            if 'ANALYSISDEC' in new_f[i]:
                dec = '-' if 'dec' not in analysis.keys() else '{:+.2f}'.format(analysis['dec'] * 180. / np.pi)
                new_f[i] = new_f[i].replace('ANALYSISDEC', dec)
            if 'ANALYSISRA' in new_f[i]:
                ra = '-' if 'ra' not in analysis.keys() else '{:.2f}'.format(analysis['ra'] * 180. / np.pi)
                new_f[i] = new_f[i].replace('ANALYSISRA', ra)
            if 'ANALYSISEXTENSION' in new_f[i]:
                ext = '-' if 'extension' not in analysis.keys() else '{:.2f}'.format(analysis['extension'] * 180. / np.pi)
                new_f[i] = new_f[i].replace('ANALYSISEXTENSION', ext)
            if 'ANALYSISNAME' in new_f[i]:
                new_f[i] = new_f[i].replace('ANALYSISNAME', analysis['name'])
            if 'ANALYSISID' in new_f[i]:
                new_f[i] = new_f[i].replace('ANALYSISID', analysis['analysisid'])    
                        
    with open('/home/apizzuto/public_html/FastResponse/webpage/output/{}.html'.format(analysis['analysisid']), 'w') as f:
        for line in new_f:
            f.write(line)

def updateFastResponseTable(analysis):
    r'''
    Push information from this analysis to overall
    tables, include new entries where applicable
    Parameters:
    -----------
    analysis: FastResponseAnalysis results from pickle file
    '''
    dec = '-' if 'dec' not in analysis.keys() else '{:+.2f}'.format(analysis['dec'] * 180. / np.pi)
    ra = '-' if 'ra' not in analysis.keys() else '{:.2f}'.format(analysis['ra'] * 180. / np.pi)
    ext = '-' if 'extension' not in analysis.keys() else '{:.2f}'.format(analysis['extension'] * 180. / np.pi)
    tag = '''
    <tr>
      <td><a href=./output/{}.html>{}</a></td>
      <td>{:.2f}</td>
      <td>{:.2e}</td>
      <td>{}</td>
      <td>{}</td>
      <td>{:.3f}</td>
    </tr>
    '''.format(analysis['analysisid'], analysis['name'], analysis['start'],
                (analysis['stop'] - analysis['start']) * 86400., 
                    ra, dec, analysis['p'])
    with open("/home/apizzuto/public_html/FastResponse/webpage/index.html", "r") as f:    
        lines = f.readlines()
    ind = None
    for i in range(len(lines)):
        if '</div></div><hr/>' in lines[i]:
            ind = i
    lines[ind-1:ind-1] = [t + '\n' for t in tag.split('\n')]
    with open("/home/apizzuto/public_html/FastResponse/webpage/index.html", 'w') as f:
        for line in lines:
            if line == '\n':
                continue
            else:
                f.write(line)

def updateFastResponsePlots():
    r'''
    Update overview plots of all analyses (timing, 
    p-value distribution, etc.)
    '''
    df = pd.read_pickle('/data/user/apizzuto/fast_response_skylab/results_dataframe.pkl')
    p_x_vals = np.logspace(-3,0.,15)
    plt.figure(figsize = (10,6), dpi=300)
    plt.hist(df['Pre-trial p_val'], weights = np.ones(len(df)) / len(df), bins = p_x_vals)
    plt.step(p_x_vals[1:], np.diff(p_x_vals), label = 'Uniform p-value expectation', lw = 3.)
    plt.xscale('log')
    plt.yscale('log')
    plt.gca().invert_xaxis()
    plt.grid(which = 'both', alpha = 0.2)
    plt.xlim(1.1e0,1e-3)
    plt.ylim(1e-2, 1e0)
    plt.xlabel('p-value', fontsize = 18)
    plt.ylabel('Fraction of Analyses', fontsize = 18)
    plt.tick_params(labelsize = 18)
    plt.legend(loc = 1, fontsize = 18)
    today = datetime.date.today().strftime("%B %d, %Y")
    plt.title("{} Fast Response Analyses as of {}".format(len(df), today), fontsize = 20)          
    #plt.text(7e-3, 5e-2, "IceCube\nPreliminary", fontsize = 20, color = 'r')
    plt.ylim(6e-3, 1e0)
    plt.savefig('/home/apizzuto/public_html/FastResponse/webpage/output/pvalue_distribution_liveupdate.png', dpi=200, bbox_inches='tight')

def write_alert_gcn(alert_results):
    r'''Read in template GCN circular file, fill in appropriate
    details, write new text file'''
    new_f = []
    high_sig = False
    for window in alert_results.keys():
        if alert_results[window]['p'] < 0.01:
            high_sig = True
    fname = 'circular_templates/internal_followup.txt'
    analysis_1000 = alert_results[1000.]; analysis_2day = alert_results[172800.]
    alert_id = analysis_1000['name'][:-10]
    if 'coincident_events' not in analysis_1000.keys():
        analysis_1000['coincident_events'] = []
    ev_is_are = 'event is' if len(analysis_1000['coincident_events']) == 1 else 'events are'
    if not high_sig:
        if len(analysis_1000['coincident_events']) == 0:
            coinc_and_p = ''
        elif len(analysis_1000['coincident_events']) == 1:
            coinc_and_p = 'We find that this additional event is well described by atmospheric\n' \
                + 'background expectations, with a p-value of {:.2f}. '.format(analysis_1000['p'])
        else:
            coinc_and_p = 'We find that these additional {} events are well described by atmospheric\n' \
                + 'background expectations, with a p-value of {:.2f}. '.format(len(analysis_1000['coincident_events']), analysis_1000['p'])
        #coinc_and_p = coinc_and_p + '\nWe accordingly derive an upper limit\n' \
        #        + 'on the time-integrated muon-neutrino flux at 1 TeV of E^2 dN/dE = {:.1e} TeV cm^-2 at 90% CL,\n'.format(analysis_1000['upper_limit']) \
        #        + 'under the assumption of an E^-2.5 power law. '
        long_p_and_lim = 'In this case, we report a p-value of {:.2f},'.format(analysis_2day['p']) \
                + ' consistent with no significant \nexcess of track events.' # , and a corresponding upper limit on the' \
                # + 'time-integrated muon-neutrino flux at 1 TeV assuming an E^-2.5 spectrum (E^2 dN/dE) of' \
                # + '{:.1e} TeV cm^-2 at the 90% CL.'.format(analysis_2day['upper_limit'])
    else:
        coinc_and_p = 'We accordingly derive a p-value of {:.3f}.'.format(analysis_1000['p'])
        if analysis_1000['p'] < 0.01:
            coinc_and_p = coinc_and_p + ' Due to the coincidences identified in this search, ' \
                + 'we strongly encourage followup observations.'
        else:
            pass
            #coinc_and_p = coinc_and_p + '\nWe accordingly derive an upper limit\n' \
            #    + 'on the time-integrated muon-neutrino flux at 1 TeV of E^2 dN/dE = {:.1e} TeV cm^-2 at 90% CL,\n'.format(analysis_1000['upper_limit']) \
            #    + 'under the assumption of an E^-2.5 power law. '
        if analysis_2day['p'] < 0.01:
            long_p_and_lim = 'In this case, we report a p-value of {:.3f}.'.format(analysis_2day['p']) \
                + ' Due to the coincidences identified in this search, we strongly encourage followup observations.'
        else:
            long_p_and_lim = 'In this case, we report a p-value of {:.2f},'.format(analysis_2day['p']) \
                + ' consistent with no significant \nexcess of track events. ' #, and a corresponding upper limit on the' \
                #+ 'time-integrated muon-neutrino flux at 1 TeV assuming an E^-2.5 spectrum (E^2 dN/dE) of' \
                #+ '{:.1e} TeV cm^-2 at the 90% CL.'.format(analysis_2day['upper_limit'])

    #if 'gcn_num' not in analysis_1000.keys():
    #    analysis_1000['gcn_num'] = 0

    keypairs = [('alert_name', alert_id), ('gcn_number', analysis_1000['gcn_num']), 
                ('start_utc', Time(analysis_1000['start'], format='mjd').iso), 
                ('stop_utc', Time(analysis_1000['stop'], format='mjd').iso), 
                ('long_start_utc', Time(analysis_2day['start'], format='mjd').iso), 
                ('long_stop_utc', Time(analysis_2day['stop'], format='mjd').iso),
                ('n_events', len(analysis_1000['coincident_events'])), 
                ('events_is_are', ev_is_are),
                ('coincident_events_and_p_str', coinc_and_p),
                ('long_p_and_lim', long_p_and_lim),
                ('low_sens', analysis_1000['sens_range'][0]), ('high_sens', analysis_1000['sens_range'][1]), 
                ('long_low_sens', analysis_2day['sens_range'][0]), ('long_high_sens', analysis_2day['sens_range'][1]),
                ('low_en', np.min([analysis_1000['low_en'], analysis_2day['low_en']])), 
                ('high_en', np.max([analysis_1000['high_en'], analysis_2day['high_en']]))]

    with open(fname, 'r') as f:
        for line in f.readlines():
            for k, r in keypairs:
                if k in line:
                    if type(r) == str:
                        form_r = r
                    elif type(r) == int:
                        form_r = '{}'.format(r)
                    elif k in ['low_en', 'high_en']:
                        form_r = '{:.0e}'.format(r)
                    elif 'sens' in k:
                        form_r = '{:.1e}'.format(r)
                    else:
                        form_r = '{}'.format(r)
                    line = line.replace('<'+k+'>', form_r)
            new_f.append(line)
    
    with open('/home/apizzuto/public_html/FastResponse/alert_circulars/{}_circular.txt'.format(alert_id), 'w') as f:
        for line in new_f:
            f.write(line)

def erfunc(x, a, b):
    x = np.array(x)
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

def binomial_error(p, number):
    errs = np.sqrt(p*(1.-p) / number)
    ntrig = p * number
    bound_case_pass = (ntrig + (1./3.)) / (number + (2./3.))
    bound_case_sigma = np.sqrt(bound_case_pass*(1. - bound_case_pass) / (number + 2))
    errs = np.maximum(errs, bound_case_sigma)
    return errs

def sync_to_roc():
    #subprocess.Popen('rsync -a /home/apizzuto/public_html/FastResponse/webpage/ apizzuto@roc.icecube.wisc.edu:/mnt/roc/www/internal/fast_response')
    env = dict(os.environ)
    subprocess.call(['rsync','-a','/home/apizzuto/public_html/FastResponse/webpage/',
                        'apizzuto@roc.icecube.wisc.edu:/mnt/roc/www/internal/fast_response'],
                        env = env
                       )

def format_dec_str(dec):
    if dec < 0:
        return str(dec)
    else:
        return '+'+str(dec)

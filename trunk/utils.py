import numpy as np
import scipy as sp
from scipy.optimize       import curve_fit
from scipy.stats          import chi2
import pandas as pd
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

def updateDataFrame(analysis):
    r'''
    Read in official Fast Response Dataframe,
    add these results, save'''
    df = pd.read_pickle('/data/user/apizzuto/fast_response_skylab/results_dataframe.pkl')
    evid = None if 'skipped' not in analysis.keys() else str(analysis['skipped']['run']) + ':' + str(analysis['skipped']['event'])
    dec = np.nan if 'dec' not in analysis.keys() else analysis['dec'] * 180. / np.pi
    ra = np.nan if 'ra' not in analysis.keys() else analysis['ra'] * 180. / np.pi
    ext = 0.0 if 'extension' not in analysis.keys() else results['extension'] * 180. / np.pi
    new_list = [ra, dec, analysis['p'],
            analysis['ns'], pd.Timestamp(Time(analysis['start'], format='mjd').iso), 
            pd.Timedelta(analysis['stop'] - analysis['start'], unit='day'),
            ext, None, analysis['ts'], evid]
    if analysis['name'] in df.index:
        num = np.count_nonzero(df.index == analysis['name'])
        analysis['name'] += '_{}'.format(num)
    df.loc[analysis['name']] = new_list
    #df.to_pickle('/data/user/apizzuto/fast_response_skylab/results_dataframe.pkl')    

def createFastResponsePage(analysis):
    r'''
    Create analysis specific page
    Parameters:
    -----------
    analysis: FastResponseAnalysis results from pickle file
    '''
    new_f = []
    keypairs = [('ANALYSISTS', 'ts'), ('ANALYSISNS', 'ns'), ('ANALYSISP', 'p')]
    with open('html/analysis_base.html', 'r') as f:
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
                ext = '-' if 'extension' not in analysis.keys() else '{:.2f}'.format(results['extension'] * 180. / np.pi)
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
    ext = '-' if 'extension' not in analysis.keys() else '{:.2f}'.format(results['extension'] * 180. / np.pi)
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
    plt.step(p_x_vals[1:], np.diff(p_x_vals), label = 'Continuous p-value expectation', lw = 3.)
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
    plt.savefig('/home/apizzuto/public_html/FastResponse/pvalue_distribution_liveupdate.png', dpi=200, bbox_inches='tight')

def write_alert_circular(analysis):
    r'''Read in template GCN circular file, fill in appropriate
    details, write new text file'''
    #keypairs = [('ANALYSISTS', 'ts'), ('ANALYSISNS', 'ns'), ('ANALYSISP', 'p')]
    if analysis['p'] < 0.01:
        fname = 'circular_templates/internal_high_significance.txt'
    else:
        fname = 'circular_templates/internal_low_significance.txt'
    keypairs = [(), (), (), (), (), (), (), (), ()]
    with open(fname, 'r') as f:
        for line in f.readlines():
            for k, r in keypairs:
                if k in line:
                    line = line.replace(k, '{:.3f}'.format(analysis[r]))
            new_f.append(line)

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

def format_dec_str(dec):
    if dec < 0:
        return str(dec)
    else:
        return '+'+str(dec)

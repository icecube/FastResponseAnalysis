import numpy as np
import scipy as sp
from scipy.optimize       import curve_fit
from scipy.stats          import chi2
import pandas as pd
import subprocess
import os, pwd
from astropy.time import Time
import datetime
import matplotlib as mpl
mpl.use('agg')
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import fast_response
import dateutil.parser
from astropy.time import Time

username = pwd.getpwuid(os.getuid())[0]
if username == 'realtime': username='jthwaites'

############################# Plotting Parameters #############################
mpl.rcParams['text.usetex'] = True
try:
    mpl.rcParams['text.latex.unicode'] = True
except:
    # new mpl doesn't like this rcParam
    pass
mpl.rcParams['mathtext.rm'] = 'Times New Roman'
mpl.rcParams['mathtext.it'] = 'Times New Roman:italic'
mpl.rcParams['mathtext.bf'] = 'Times New Roman:bold'

mpl.rc('font', family='serif', size=12)
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['ytick.major.size'] = 5
############################# Plotting Parameters #############################


def updateFastResponseWeb(analysis, gw=False):
    r'''
    Create analysis specific page, and update
    plots with information from all analyses
    Parameters:
    -----------
    analysis: FastResponseAnalysis instance
    '''
    updateDataFrame(analysis, gw=gw)
    createFastResponsePage(analysis, gw=gw)
    if gw:
        updateGWTable(analysis)
    else:
        updateFastResponseTable(analysis)
    updateFastResponsePlots(gw=gw)
    # sync_to_roc()

def updateDataFrame(analysis, gw=False, make_df=False):
    r'''
    Read in official Fast Response Dataframe,
    add these results, save
    gw option specifies to read in GW DataFrame
    make_df allows to make a new dataframe with specific keys in same format
    '''
    base_path = os.path.dirname(fast_response.__file__)
    if make_df:
        if gw:
            print('Making new dataframe')
            df=pd.DataFrame({'Pre-trial p_val':[], 
                    'Beginning Observation Time (UTC)':[],'Duration':[],
                    'Fit ns':[], 'TS':[], 'fit gamma':[], 
                    'Best-fit RA':[],'Best-fit Dec':[],
                    'Energy range':[], 'Sensitivity Range':[]})
            df.index.name='Source Name'
    else: 
        if gw: 
            df = pd.read_pickle(f'{base_path}/results_dataframe_gw.pkl')
        else: 
            df = pd.read_pickle(f'{base_path}/results_dataframe.pkl')
    #df = pd.read_pickle(f'{base_path}/results_dataframe.pkl')

    evid = None if 'skipped' not in analysis.keys() else str(analysis['skipped']['run']) + ':' + str(analysis['skipped']['event'])
    dec = np.nan if 'dec' not in analysis.keys() else analysis['dec'] * 180. / np.pi
    ra = np.nan if 'ra' not in analysis.keys() else analysis['ra'] * 180. / np.pi
    ext = 0.0 if 'extension' not in analysis.keys() else analysis['extension'] * 180. / np.pi
    upper_lim = np.nan if 'upper_limit' not in analysis.keys() else analysis['upper_limit']
    #note that this upper_lim will be np.nan for any skymap followup, because then we send out a sens_range
    
    #low_en = np.nan if 'low_en' not in analysis.keys() else analysis['low_en']
    #high_en = np.nan if 'high_en' not in analysis.keys() else analysis['high_en']

    if gw:
        new_list = [analysis['p'], pd.Timestamp(Time(analysis['start'], format='mjd').iso),
            pd.Timedelta(analysis['stop'] - analysis['start'], unit='day'),
            analysis['ns'], analysis['ts'], analysis['gamma'], analysis['fit_ra'], 
            analysis['fit_dec'], analysis['energy_range'], analysis['sens_range']]
    else: 
        new_list = [ra, dec, analysis['p'],
                analysis['ns'], pd.Timestamp(Time(analysis['start'], format='mjd').iso), 
                pd.Timedelta(analysis['stop'] - analysis['start'], unit='day'),
                ext, None, analysis['ts'], evid, upper_lim, analysis['energy_range']]
    if analysis['name'] in df.index:
        num = np.count_nonzero(df.index == analysis['name'])
        analysis['name'] += '_{}'.format(num)
    df.loc[analysis['name']] = new_list
     
    if gw: 
        df.to_pickle(f'{base_path}/results_dataframe_gw.pkl')
    else:
        df.to_pickle(f'{base_path}/results_dataframe.pkl')

def createFastResponsePage(analysis, gw=False):
    r'''
    Create analysis specific page
    Parameters:
    -----------
    analysis: FastResponseAnalysis results from pickle file
    gw: specifier if analysis is for a GW followup
    '''
    new_f = []
    keypairs = [('ANALYSISTS', 'ts'), ('ANALYSISNS', 'ns'), ('ANALYSISP', 'p')]
    base_path = os.path.dirname(fast_response.__file__)
    html_base = f'{base_path}/../html/'
    
    with open(f'{html_base}analysis_base.html', 'r') as f:
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
                new_f[i] = new_f[i].replace('ANALYSISDURATION', str(round(dur,2)))
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
            if 'ANALYZER' in new_f[i]:
                new_f[i] = new_f[i].replace('ANALYZER', username)
            if gw:
                if 'webpage' in new_f[i]:
                    new_f[i] = new_f[i].replace('webpage', 'gw-webpage')
                if '<tr><td>p-value:</td><td>' in new_f[i]:
                    ind=i

    if gw: 
        gracedb_link = 'https://gracedb.ligo.org/superevents/{}/view'.format(analysis['name'].split('-')[0])
        new_f[ind+1:ind+1] = '    <tr><td>GraceDB link:</td><td><a href={}>link</a></td></tr>\n'.format(gracedb_link)
        webpage_path='/home/{}/public_html/FastResponse/gw-webpage/output/{}.html'.format(username, analysis['analysisid'])
    else: 
        webpage_path='/home/{}/public_html/FastResponse/webpage/output/{}.html'.format(username, analysis['analysisid'])
    with open(webpage_path, 'w') as f:
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
    
    with open(f"/home/{username}/public_html/FastResponse/webpage/index.html", "r") as f:    
        lines = f.readlines()
    ind = None
    for i in range(len(lines)):
        if '</div></div><hr/>' in lines[i]:
            ind = i
    lines[ind-1:ind-1] = [t + '\n' for t in tag.split('\n')]
    with open(f"/home/{username}/public_html/FastResponse/webpage/index.html", 'w') as f:
        for line in lines:
            if line == '\n':
                continue
            else:
                f.write(line)

def updateGWTable(analysis):
    r'''
    Push information from this analysis to summary tables 
    on GW followup webpage
    Parameters:
    -----------
    analysis: GW FastResponseAnalysis results from pickle file
    '''
    dec = '-' if analysis['ts']<=0. else '{:+.2f}'.format(analysis['fit_dec'] * 180. / np.pi)
    ra = '-' if analysis['ts']<=0. else '{:.2f}'.format(analysis['fit_ra'] * 180. / np.pi)
    ts=analysis['ts']
    #ts= max(0.,analysis['ts'])

    tag = '''
    <tr>
      <td><a href=./output/{}.html>{}</a></td>
      <td>{:.2f}</td>
      <td>{}</td>
      <td>{}</td>
      <td>{:.3f}</td>
      <td>{:.3f}</td>
      <td>{:.2f}</td>
      <td>{:.3f}</td>
    </tr>
    '''.format(analysis['analysisid'], analysis['name'], 
                analysis['start'], ra, dec, ts,
                analysis['ns'], analysis['gamma'], 
                analysis['p'])
    
    with open(f"/home/{username}/public_html/FastResponse/gw-webpage/index.html", "r") as f:    
        lines = f.readlines()
    ind = None
    for i in range(len(lines)):
        if '</div></div><hr/>' in lines[i]:
            ind = i
    lines[ind-1:ind-1] = [t + '\n' for t in tag.split('\n')]
    with open(f"/home/{username}/public_html/FastResponse/gw-webpage/index.html", 'w') as f:
        for line in lines:
            if line == '\n':
                continue
            else:
                f.write(line)

def updateFastResponsePlots(gw=False):
    r'''
    Update overview plots of all analyses (timing, 
    p-value distribution, etc.)
    '''
    base_path = os.path.dirname(fast_response.__file__)

    if gw: 
        df = pd.read_pickle(f'{base_path}/results_dataframe_gw.pkl')
    else: 
        df = pd.read_pickle(f'{base_path}/results_dataframe.pkl')
    #df = pd.read_pickle(f'{base_path}/results_dataframe.pkl')
    
    p_x_vals = np.logspace(-3,0.,15)
    plt.figure(figsize = (10,6), dpi=300)
    plt.hist(df['Pre-trial p_val'], weights = np.ones(len(df)) / len(df), bins = p_x_vals)
    plt.step(p_x_vals[1:], np.diff(p_x_vals), label = 'Uniform p-value distribution', lw = 3.)
    plt.xscale('log')
    plt.yscale('log')
    plt.gca().invert_xaxis()
    plt.grid(which = 'both', alpha = 0.2)
    plt.xlim(1.1e0,1e-3)
    plt.ylim(3e-3, 1e0)
    plt.xlabel('p-value', fontsize = 18)
    plt.ylabel('Fraction of Analyses', fontsize = 18)
    plt.tick_params(labelsize = 18)
    plt.legend(loc = 1, fontsize = 18)
    today = datetime.date.today().strftime("%B %d, %Y")

    if gw: 
        pval_dist_path=f'/home/{username}/public_html/FastResponse/gw-webpage/output/pvalue_distribution_liveupdate.png'
        plt.title("{} GW Followups as of {}".format(len(df), today), fontsize = 20) 
    else:
        pval_dist_path=f'/home/{username}/public_html/FastResponse/webpage/output/pvalue_distribution_liveupdate.png'
        plt.title("{} Fast Response Analyses as of {}".format(len(df), today), fontsize = 20)          
    #plt.text(7e-3, 5e-2, "IceCube\nPreliminary", fontsize = 20, color = 'r')
    plt.ylim(3e-3, 1e0)
    
    plt.savefig(pval_dist_path, dpi=200, bbox_inches='tight')
    #plt.savefig(f'/home/{username}/public_html/FastResponse/webpage/output/pvalue_distribution_liveupdate.png', dpi=200, bbox_inches='tight')

def updateGW_public(analysis, circular = None):
    r'''
    Push information from this analysis to summary tables 
    on GW PUBLIC followup webpage
    Parameters:
    -----------
    analysis: dictionary of results from UML and LLAMA analyses
    circular: GCN circular number, to link to circular rather than notice
    '''

    #duration = (Time(analysis['observation_stop'], format='iso').mjd - Time(analysis['observation_start'],format = 'iso').mjd)*86400.
    duration = analysis['observation_livetime']
    if 'analysisid' not in analysis.keys():
        start_date = Time(dateutil.parser.parse(analysis['observation_start'])).datetime
        start_str = f'{start_date.year:02d}_' \
            + f'{start_date.month:02d}_{start_date.day:02d}'
        analysis['analysisid'] = start_str+'_'+analysis['reference']['gcn.notices.LVK.alert']

    analysis['name'] = analysis['reference']['gcn.notices.LVK.alert']
    if 'test' in analysis['analysisid']:
        analysis['name'] = analysis['name']+'_test'

    subprocess.call(['cp', 
                     '/home/followup/lvk_followup_output/{}_collected_results.json'.format(analysis['name']),
                     f'/home/{username}/public_html/public_FRA/gw-webpage/output/'])
    
    if analysis['n_events_coincident'] == 0:
        extra_info = '-'
    else:
        outdir = os.path.join(os.environ.get('FAST_RESPONSE_OUTPUT'),analysis['analysisid'])
        try:
            if not analysis['subthreshold'] and analysis['pval_generic']<0.1:
                subprocess.call(['cp', '{}/{}unblinded_skymap_zoom.png'.format(outdir,analysis['analysisid']),
                             f'/home/{username}/public_html/public_FRA/gw-webpage/output/'])
        except:
            print('Failed to copy skymap and GCN Notice')
        createGWEventPage(analysis)
        extra_info = '<a href=./output/{}.html>here</a>'.format(analysis['analysisid'])

    if 'skymap_path' not in analysis.keys():
        analysis['skymap_path'] = ''
    if 'subthreshold' not in analysis.keys():
        analysis['subthreshold']=False
    if circular is not None:
        link = '<a href=https://gcn.nasa.gov/circulars/{}>GCN Circular</a>'.format(circular)
    else:
        link = '<a href=./output/{}_collected_results.json>GCN Notice</a>'.format(analysis['name'])

    if analysis['subthreshold']:
        row_start = '<tr bgcolor="#ddd" name="subthr">'
    else:
        row_start = '<tr name="sig">'
    tag = '''
    {}
        <td>{}</td>
        <td>{}</td>
        <td><a href=https://gcn.gsfc.nasa.gov/notices_l/{}.lvc>LVK GCN</a></td>
        <td>{}</td>
        <td>{:.2e}</td>
        <td>{}</td>
        <td>[{:.3f}, {:.3f}]</td>        
        <td>{}</td>
    </tr>
    '''.format(row_start,analysis['name'], analysis['trigger_time'][:-1].replace('T',' '),
               analysis['name'].split('-')[0], link, duration, 
               str(analysis['n_events_coincident']), 
               analysis["neutrino_flux_sensitivity_range"]['flux_sensitivity'][0],
               analysis["neutrino_flux_sensitivity_range"]['flux_sensitivity'][1], extra_info)
    
    with open(f"/home/{username}/public_html/public_FRA/gw-webpage/index.html", "r") as f:    
        lines = f.readlines()
    ind = None
    for i in range(len(lines)):
        if '</thead>' in lines[i]:
            ind = i
    lines[ind+1:ind+1] = [t + '\n' for t in tag.split('\n')]
    with open(f"/home/{username}/public_html/public_FRA/gw-webpage/index.html", 'w') as f:
        for line in lines:
            if line == '\n':
                continue
            else:
                f.write(line)

def createGWEventPage(analysis):
    r'''
    Create PUBLIC webpage with neutrino information, if sent via GCN
    Parameters:
    -----------
    analysis: dictionary of results from UML and LLAMA analyses
    '''
    new_f = []
    base_path = os.path.dirname(fast_response.__file__)
    html_base = f'{base_path}/../html/'

    with open(f'{html_base}gw_pub_templ_base.html', 'r') as f:
        for line in f.readlines():
            new_f.append(line)
        ind = None
        for i in range(len(new_f)):
            if 'ANALYSISCREATIONTIME' in new_f[i]:
                new_f[i] = new_f[i].replace('ANALYSISCREATIONTIME', str(datetime.datetime.utcnow())[:-7])
            if 'ANALYSISSTART' in new_f[i]:
                new_f[i] = new_f[i].replace('ANALYSISSTART', analysis['observation_start'])
            if 'ANALYSISSTOP' in new_f[i]:
                new_f[i] = new_f[i].replace('ANALYSISSTOP', analysis['observation_stop'])
            if 'ANALYSISNAME' in new_f[i]:
                new_f[i] = new_f[i].replace('ANALYSISNAME', analysis['name'])
            if 'ANALYSISID' in new_f[i]:
                new_f[i] = new_f[i].replace('ANALYSISID', analysis['analysisid'])  
            if 'NEVENTS' in new_f[i]:
                new_f[i] = new_f[i].replace('NEVENTS', str(analysis['n_events_coincident']))
            if 'UMLPVAL' in new_f[i]:
                try:
                    new_f[i] = new_f[i].replace('UMLPVAL', '{}'.format(analysis['pval_generic']))
                except:
                    new_f[i] = new_f[i].replace('UMLPVAL', 'nan')
            #if 'UMLSIG' in new_f[i]:
            #    new_f[i] = new_f[i].replace('UMLSIG', '{:.2f}'.format(analysis['overall_sig_gen_transient']))
            if 'LLAMAPVAL' in new_f[i]:
                try: 
                    new_f[i] = new_f[i].replace('LLAMAPVAL', '{}'.format(analysis['pval_bayesian']))
                except:
                    new_f[i] = new_f[i].replace('LLAMAPVAL', 'nan')
            #if 'LLAMASIG' in new_f[i]:
            #    new_f[i] = new_f[i].replace('LLAMASIG', '{:.2f}'.format(analysis['overall_sig_bayesian']))
            if 'MOSTPROBDIR' in new_f[i]:
                if 'most_probable_direction' in analysis.keys():
                    new_f[i] = new_f[i].replace('MOSTPROBDIR', 'RA: {} deg, decl.: {} deg'.format(
                                                analysis['most_probable_direction']['ra'], 
                                                analysis['most_probable_direction']['dec']))
                else: 
                    new_f[i] = new_f[i].replace('MOSTPROBDIR', 'N/A')
            if '</div></div>' in new_f[i]:
                ind = i
    
    e=1
    tag=''
    for event in analysis['coincident_events']:
        if 'event_pval_bayesian' not in event.keys():
            event['event_pval_bayesian'] = 'nan'
        if 'event_pval_generic' not in event.keys():
            event['event_pval_generic'] = 'nan'
        tag += '''
        <h3>Event {}</h3>
        <table class="skymaps">
            <table class='event'>
                <tr><td>dt (sec)</td><td>{:.2f}</td></tr>
                <tr><td>RA (J2000) (deg.):</td><td>{:.2f}</td></tr>
                <tr><td>Dec (J2000) (deg.):</td><td>{:.2f}</td></tr>
                <tr><td>Angular Unc (deg.):</td><td>{:.2f}</td></tr>
                <tr><td>Event p-value (generic transient):</td><td>{}</td></tr>
                <tr><td>Event p-value (Bayesian):</td><td>{}</td></tr>
            </table>
        </table>
        '''.format(e, event['event_dt'], event['localization']['ra'], event['localization']['dec'], 
                   event['localization']['ra_uncertainty'], 
                   event['event_pval_generic'], event['event_pval_bayesian'])
        e+=1
    
    if not analysis['subthreshold'] and analysis['pval_generic']<0.1:
        tag2 = '''
        <div class="content">
        <a href="./{}unblinded_skymap_zoom.png"><img src="./{}unblinded_skymap_zoom.png" width="500"/></a>
        '''.format(analysis['analysisid'], analysis['analysisid'])
        new_f[ind+1:ind+1] = [t + '\n' for t in tag2.split('\n')]

    new_f[ind:ind] = [t + '\n' for t in tag.split('\n')]

    webpage_path='/home/{}/public_html/public_FRA/gw-webpage/output/{}.html'.format(username, analysis['analysisid'])
    with open(webpage_path, 'w') as f:
        for line in new_f:
            f.write(line)

def update_internal_public(analysis_1000, analysis_2d, alert_gcn=None, fra_circular=None):
    r'''
    Push information from this analysis to summary tables 
    on PUBLIC FRA followup webpage
    Parameters:
    -----------
    analysis: dictionary of results
    alert_gcn: either a single int (Circular number) or tuple (run id, event id) if linking to notice
    fra_circular: GCN circular number for FRA followup
    '''
    event_name = analysis_2d['name'].replace(' 1.7e+05 s','')
    cascade=True if 'Cascade' in event_name else False

    if type(alert_gcn) == int:
        #link to circular for alert event
        alert_link = 'href=https://gcn.nasa.gov/circulars/{}'.format(alert_gcn, event_name)
    elif type(alert_gcn) == list or type(alert_gcn)==tuple:
        if len(alert_gcn)!=2:
            print('If linking to notice, need both run id and event id!')
        if cascade:
            alert_link = 'href=https://gcn.gsfc.nasa.gov/notices_amon_icecube_cascade/{}_{}.amon'.format(alert_gcn[0],alert_gcn[1])
        else:
            alert_link = 'href=https://gcn.gsfc.nasa.gov/notices_amon_g_b/{}_{}.amon'.format(alert_gcn[0],alert_gcn[1])
    else:
        alert_link=''
        print('Link to alert GCN missing!')
    event_link = '<a {}>{}</a>'.format(alert_link, event_name)

    if fra_circular is not None:
        fra_link = '<a href=https://gcn.nasa.gov/circulars/{}>{}</a>'.format(fra_circular,fra_circular)
    else:
        fra_link=' '
    
    if '{:.1e}'.format(analysis_1000['sens_range'][0]) == '{:.1e}'.format(analysis_1000['sens_range'][1]):
        sens_range_1000 = f'{analysis_1000["sens_range"][0]:.1e}'
    else:
        sens_range_1000 = f'[{analysis_1000["sens_range"][0]:.1e}, {analysis_1000["sens_range"][1]:.1e}]'

    if '{:.1e}'.format(analysis_2d['sens_range'][0]) == '{:.1e}'.format(analysis_2d['sens_range'][1]):
        sens_range_2d = f'{analysis_2d["sens_range"][0]:.1e}'
    else:
        sens_range_2d = f'[{analysis_2d["sens_range"][0]:.1e}, {analysis_2d["sens_range"][1]:.1e}]'

    tag = '''
        <tbody>
        {}
            <td rowspan='2'>{}</td>
            <td rowspan='2'>{}</td>
            <td rowspan='2'>{}</td>
            <td rowspan='2'>{}</td>
            <td>1000 seconds</td>
            <td>{}</td>
            <td rowspan='2'>2.5</td>
            <td rowspan='2'>{}</td>
        </tr>
        {}
            <td>2 days</td>
            <td>{}</td>
        </tr>
        </tbody>
    '''.format('<tr name="casc">' if cascade else '<tr name="track">',
               event_link, fra_link, Time(analysis_2d['start']+1, format='mjd').iso,
               'Cascade' if cascade else 'Track', sens_range_1000, 
               f'[{analysis_1000["energy_range"][0]:.0e}, {analysis_1000["energy_range"][1]:.0e}]',
               '<tr name="casc">' if cascade else '<tr name="track">',
               sens_range_2d)
    
    with open(f"/home/{username}/public_html/public_FRA/internal_alerts_webpage/index.html", "r") as f:    
        lines = f.readlines()
    ind = None
    for i in range(len(lines)):
        if '</thead>' in lines[i]:
            ind = i
    lines[ind+1:ind+1] = [t + '\n' for t in tag.split('\n')]
    with open(f"/home/{username}/public_html/public_FRA/internal_alerts_webpage/index.html", 'w') as f:
        for line in lines:
            if line == '\n':
                continue
            else:
                f.write(line)

#def sync_to_roc():
#    #subprocess.Popen('rsync -a /home/apizzuto/public_html/FastResponse/webpage/ apizzuto@roc.icecube.wisc.edu:/mnt/roc/www/internal/fast_response')
#    env = dict(os.environ)
#    subprocess.call(['rsync', '-a', f'/home/{username}/public_html/FastResponse/webpage/',
#                        f'{username}@roc.icecube.wisc.edu:/mnt/roc/www/internal/fast_response'],
#                        env = env
#                       )

def format_dec_str(dec):
    if dec < 0:
        return str(dec)
    else:
        return '+'+str(dec)

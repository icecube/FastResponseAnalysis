import seaborn as sns
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
import icecube.realtime_tools.live
import matplotlib as mpl
mpl.use('agg')
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

############################# Plotting Parameters #############################
mpl.use('agg')
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

current_palette = sns.color_palette('colorblind', 10)
############################# Plotting Parameters #############################


def time_axis(ax, run_table, time_window):
    """Make an axis that is in times, formatted as ISO times

    Parameters
    -----------
    ax: matplotlib.pyplot Axis
        axis to convert to times and format
    run_table: dict-like
        Table of run information for events
    time_window: tuple astropy Time object
        (start, stop) of the time window, in mjd
    """
    for run in run_table:
        plt.axvline(Time(run['start'],format='iso',scale='utc').plot_date, ls = '--', lw = 0.75, c = 'grey')
        plt.axvline(Time(run['stop'],format='iso',scale='utc').plot_date, ls = '--', lw = 0.75, c = 'grey')
        
    plt.axvspan(time_window[0].plot_date, time_window[1].plot_date, alpha = 0.2)

    if (Time(run_table[-1]['stop'],format='iso',scale='utc').mjd - 
        Time(run_table[0]['start'],format='iso',scale='utc').mjd) > 41.:
        
        #longer followups: need less ticks to plot
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=10))
        ax.xaxis.set_major_formatter( mdates.DateFormatter("%m"))
    else:
        ax.xaxis.set_minor_locator(mdates.HourLocator())
        ax.xaxis.set_major_locator(mdates.HourLocator([0,12]))
        ax.xaxis.set_major_formatter( mdates.DateFormatter("%m-%d %H:%M"))

    ax.set_xlim(Time(run_table[0]['start'],format='iso',scale='utc').plot_date,
                Time(run_table[-1]['stop'],format='iso',scale='utc').plot_date)

def time_series(ax, run_table, time_window, t1, t2, n, 
                scale=1, ymax=0, xerr = False, **kwargs):
    """Make a plot that is a rate over time

    Parameters
    -----------
    ax: matplotlib.pyplot Axis
        axis to convert to times and format
    run_table: dict-like
        Table of run information for events
    time_window: tuple astropy Time object
        (start, stop) of the time window, in mjd
    t1: astropy Time object
        start time to use in the plot, utc
    t2: astropy Time object
        stop time to use in the plot, utc
    n: list or array-like
        value to plot on the y-axis (usually, a rate value)
    scale: float
        rescale value/unit for the y-axis (default 1, not rescaled)
    ymax: float
        maximum value to plot for the y-axis (default 0, no max)
    kwargs: optional
        passed to matplotlib.pyplot.ax.errorbar
    """
    if ax is None:
        ax=plt.gca()
    time_axis(ax, run_table, time_window)

    tdiff = t2 - t1
    tmid = tdiff/2 + t1
    
    r = n/tdiff.sec*scale
    rerr = np.sqrt(n)/tdiff.sec*scale
    m = np.isfinite(r)
    
    ax.errorbar(tmid.plot_date, r,
                    xerr=tdiff.jd/2, yerr=rerr, capsize = 4,
                    ls='', c = current_palette[1], **kwargs)
    
    Ymax = max(max(r[m])*1.5,ymax)
    ax.set_ylim(-0.5,Ymax)

def make_rate_plots(time_window, run_table, query_events, dirname, season='neutrino'):
    """ Make plots of filter rates for the time window we're interested in. 
    Rates plotted are: 

    * GFU singlet rate
    * Badness
    * Online L2
    * Muon filter
    * In Ice Simple Multiplicity

    Parameters
    -----------
    time_window: tuple astropy Time object
        (start, stop) of the time window, in mjd
    run_table: dict-like
        Table of run information for events
    dirname: str
        Directory events are saved to
    season: str
        Season to load (default neutrino)
    """
    ########## MAKE GFU RATE PLOT ##########
    fig, ax = plt.subplots(figsize = (12,4))

    gfu_run_start = Time([ r['start'] for r in run_table], format='iso',scale='utc')
    gfu_run_stop  = Time([r['stop' ] for r in run_table], format='iso',scale='utc')
    gfu_count = [r['gfu_counts'] for r in run_table]

    time_series(ax, run_table, time_window, 
                gfu_run_start, 
                gfu_run_stop, 
                gfu_count, scale = 1e3)

    plt.title('GFU Singlet Rate', fontsize = 18)
    plt.ylabel('Rate (mHz)', fontsize = 16)
    fig.autofmt_xdate()
    plt.locator_params(axis='x', nbins = 8)
    plt.grid(b = True, axis = 'y', alpha = 0.3)
    plt.savefig('{}/GFU_rate_plot.png'.format(dirname))

    try:
        badness = icecube.realtime_tools.live.get_badness(run_table[0]['start'], run_table[-1]['stop'])
    except Exception as e:
        print(e)
        print("Failed to make badness plot. Failed to load badness values (line 95, make_ontime_plots.py)")
        return

    ########## MAKE BADNESS PLOT ##########   
    fig, ax = plt.subplots(figsize = (12,4))
    time_series(ax, run_table, time_window, 
                Time(badness['start'],format='mjd',scale='utc'), 
                Time(badness['stop' ],format='mjd',scale='utc'), 
                badness['score'])

    plt.title('Badness', fontsize = 18)
    plt.ylabel('Badness (s)', fontsize = 16)
    fig.autofmt_xdate()
    plt.locator_params(axis='x', nbins = 8)
    plt.grid(b = True, axis = 'y', alpha = 0.3)
    plt.savefig('{}/badness_plot.png'.format(dirname))
     
    try:
        rates = icecube.realtime_tools.live.get_rates(run_table[0]['start'], run_table[-1]['stop'])
    except Exception as e:
        #rates=np.load('/data/user/jthwaites/FastResponseAnalysis/output/2022_2mo_followup_rates.npy')
        print(e)
        print("Failed to load rates. In case of timeout: save rate information and reload manually")
        print("Check lines 115-119 of make_ontime_plots.py")
        return
        
    ########## MAKE RATES PLOTS ##########
    recstart= Time(rates['rec_start'],format='mjd',scale='utc')
    recstop = Time(rates['rec_stop' ],format='mjd',scale='utc')

    if season == 'neutrino16':
        online_str = 'OnlineL2Filter_16' 
        muon_str='MuonFilter_13'
    elif time_window[1].mjd > 60276.86875:
        online_str='OnlineL2Filter_23'
        muon_str='MuonFilter_23'
    else:
        online_str = 'OnlineL2Filter_17'
        muon_str='MuonFilter_13'

    for rate, name, unit in [('IN_ICE_SIMPLE_MULTIPLICITY', 'In-Ice-Simple-Multiplicity', '(kHz)'), 
                (muon_str, 'Muon Filter', '(Hz)'), 
                (online_str, 'Online L2 Filter', '(Hz)')]:
        fig, ax = plt.subplots(figsize = (12,4))
        time_series(ax, run_table, time_window, 
                    recstart,  
                    recstop, 
                    rates[rate])
            
        plt.title(name, fontsize = 18)
        plt.ylabel('{} {}'.format(name, unit), fontsize = 16)
        fig.autofmt_xdate()
        plt.locator_params(axis='x', nbins = 8)
        plt.grid(b = True, axis = 'y', alpha = 0.3)
        plt.savefig('{}/{}_plot.png'.format(dirname, rate))
        

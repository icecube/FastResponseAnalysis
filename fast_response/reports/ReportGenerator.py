'''Class for generating reports for Fast Reponse and GW
   followup analyses. A lot of code has been taken from Fast Reponse
   Analysis. 

   Author: Kevin Meagher, Raamis Hussain & Alex Pizzuto
   Date:   Mar 27, 2019
   Modified: April, 2019
'''

import numpy as np
import pandas as pd
import os, json, datetime
import logging as log

use_urllib2 = True
try:
    import urllib2, urllib
except:
    use_urllib2 = False
    import urllib.request, urllib.parse, urllib.error

import icecube.realtime_tools.live
import subprocess
from os.path import expanduser

from icecube             import icetray
import skylab
from skylab.datasets import Datasets
from astropy.time        import Time, TimeDelta
from fast_response.make_ontime_plots   import make_rate_plots
import fast_response

log.basicConfig(level=log.ERROR)
mpl_logger = log.getLogger('matplotlib') 
mpl_logger.setLevel(log.ERROR) 
 
class ReportGenerator(object):
    r"""
    Generate report from a FastResponseAnalysis object
    """
    _figure_list = []

    def __init__(self, analysis):
        self.analysis = analysis
        if self.analysis.skymap is None:
            self.source_type = 'PS'
        else:
            self.source_type = 'skymap'

        source = {}
        source['trigger_iso'] = Time(analysis.centertime, format='mjd').iso
        source['start_iso'] = Time(analysis.start, format='mjd').iso
        source['stop_iso'] = Time(analysis.stop, format='mjd').iso
        source['realtime'] = (analysis.stop - analysis.start) * 86400.
        source['start_mjd'] = analysis.start
        source['stop_mjd'] = analysis.stop
        source['trigger_mjd'] = analysis.centertime
        self.source = source

        if self.source["stop_mjd"] > 58309.74:
            self.stream = 'neutrino'
        elif self.source["stop_mjd"] > 57891.17:
            self.stream = 'neutrino17'
        else:
            self.stream = 'neutrino16'

        self.analysisid = analysis.analysisid
        self.dirname = analysis.analysispath
        self.time_window=(
            Time(self.source['start_mjd'], format='mjd'), 
            Time(self.source['stop_mjd'], format='mjd'))
        self.trigger = Time(self.source['trigger_mjd'], format='mjd')

    def get_report_source(self):
        """ Find base path of the general report TeX document
        Looks for the directory location fast_response, returns 
        fast_response/latex/report_general.tex
        """
        base = os.path.dirname(fast_response.__file__)
        return os.path.join(base, 'latex/report_general.tex')

    def write_table(self, file, name, header, table, prefix=""):
        """
        Write a table in latex format for the generated reports

        Parameters
        -----------
        file: I/O file object
            Opened file that is being written to
        name: str
            Name of the table
        header: list
            Headers for each column (can be empty if none)
        table: array of tuples
            set of values to put in the table. tuples give individual rows
        prefix: str
            prefix to latex command, if needed (default "")
        """

        cols = max([len(r) for r in table if r is not None]+[len(header)])

        x=''.join(' & '.join(str(x).strip() for x in row)+r'\\'+'\n' for row in table if row is not None)

        file.write(
            prefix+
            r"\newcommand{"+"\\"+name+"}{"+"\n"+
            r"\begin{longtable}{" + 'l' * cols + r"}"+'\n'+
            r"\hline"+'\n'+
            ' & '.join(str(x).strip() for x in header)+
            (r"\\ \hline"+'\n' if header else "") +
            ''.join(' & '.join(str(x).strip() for x in row)+r'\\'+'\n' if row is not None else r"\hline"+'\n' for row in table)+
            ((r"\hline"+'\n') if table and table[-1] is not None else "") +
            r"\end{longtable}"+
            "}\n")

    def query_db_runs(self, time_window):
        """ Queries IceCube Live for the ontime window, 
        plus 8 hours on either side (or as much as possible, if in realtime)

        Parameters
        -----------
        time_window: tuple of astropy.time Time objects
            (start, stop) for the analysis, in mjd
        
        Returns
        ----------
        dict: 
            table of runs loaded, and status of the runs
        """
        
        run_url = 'https://live.icecube.wisc.edu/run_info/'
        with open('/home/jthwaites/private/tokens/auth.txt') as f:
            user = f.readline().rstrip('\n')
            pasw = f.readline().rstrip('\n')
        run_query = {'user':user, 'pass':pasw,
                     'start':(Time(time_window[0],precision=0)
                              - TimeDelta(8*3600,format='sec')).iso,
                     'stop': (Time(time_window[1],precision=0)
                              + TimeDelta(8*3600,format='sec')).iso}
        now = Time(datetime.datetime.now()+datetime.timedelta(hours=6),scale='utc',precision=0)
        if use_urllib2:
            run_table = json.loads(urllib2.urlopen(
                urllib2.Request(run_url, urllib.urlencode(run_query)),
                timeout=500).read())
        else:
            req_data = urllib.parse.urlencode(run_query).encode("utf-8")
            req = urllib.request.Request(run_url)
            with urllib.request.urlopen(req, data=req_data, timeout=500) as fi:
                run_table = json.load(fi)

        #add duration to run table
        for run in run_table:
            runstart = Time(run['start'],format='iso',scale='utc')
            if run['stop'] is None:
                runstop  = now
                run['stop']=runstop.iso
            else:
                runstop  = Time(run['stop'],scale='utc',precision=0)

            dt=int((runstop-runstart).sec)
            hours = int(dt)//3600
            minutes = int(dt)//60-hours*60
            seconds = dt%60

            livetime = (min(time_window[1],runstop)-
                        max(time_window[0],runstart)).sec
            run['livetime'] = livetime if livetime > 0 else 0
            if (run['livetime'] > 86400. * 2) or (hours > 100):
                run['livetime'] = 0
                run['stop'] = run['start']
                dt=0
                hours = 0
                minutes = 0
                seconds = 0

            run['duration_sec'] = dt
            run['duration'] = "{:01d}:{:02d}:{:02d}".format(hours,minutes,seconds)
            run['OK'] = ("OK" if run['status'] in ['SUCCESS']
                         and run['lightmode'] in ['dark']
                         and run['filter_mode'] in ['PhysicsFiltering']
                         and run['run_mode'] in ['PhysicsTrig']
                         else "NotOK"
                         )

        return run_table

    @staticmethod
    def ontime_table(query_dict):
        """Makes a DataFrame of events, queried from IceCube Live

        Parameters
        -----------
        query_dict: dict
            Queried events in a particular time window
        
        Returns
        ----------
        pandas DataFrame:
            Event parameters for the given events
        """
        newdict=[]
        for event in query_dict:
            newevent = event['value']['data']
            for key,val in list(newevent['reco']['splinempe'].items()): 
                newevent['splinempe_'+key]=val
            if '+' in newevent['eventtime']: 
                newevent['eventtime'] = newevent['eventtime'].split('+')[0]
            if Time(newevent['eventtime'],scale='utc',format='iso') >= Time("2018-07-10 17:52:03.34", format='iso',scale='utc'):
                newevent['muex'] = newevent['reco']['energy']['mpe_muex']
            del newevent['reco']
            newdict.append(newevent)
        events = pd.DataFrame(newdict)

        if len(events):
            t=Time(list(events['eventtime']),scale='utc',format='iso')
            events['t']=t.datetime
            events['mjd'] = t.mjd
        return events

    def make_coinc_events_table(self, f):
        """Make a table of ontime events, in LaTeX format

        Parameters
        -----------
        f: I/O file object
            Opened file that is being written to
        """
        event_table = []
        if self.analysis.coincident_events is not None and self.analysis.coincident_events != []:
            if self.analysis.skymap is None:
                for event in self.analysis.coincident_events:
                    event_table+=[
                        ("Run:Event",'{}:{}'.format(event['run'], event['event'])),
                        ("Time","{}".format(
                            event['time'])),
                        (r'$\alpha$, $\delta$',"{:3.2f}\degree, {:+3.2f}\degree"
                            .format(np.rad2deg(event['ra']), np.rad2deg(event['dec']))),
                        ("Angular Uncertainty (90\%)","{:3.2f}\degree"
                            .format(np.rad2deg(event["sigma"]*2.145966))),
                        ("Distance from Source", "{:2.2f}\degree".format(event['delta_psi']*180. / np.pi)),
                        ("Reconstructed Energy (GeV)","{:2.2e}"
                            .format(10**event['logE'])),
                        ("Spatial Weight", "{:.2f}".format(event['spatial_w'])),
                        ("Energy Weight", "{:.2f}".format(event['energy_w'])),
                        None,
                    ]

                self.write_table(f, "event", [], event_table)
            else:
                if 'pvalue' in self.analysis.coincident_events[0].keys():
                    with_p = True
                    if len(self.analysis.coincident_events)<100:
                        event_table+=[
                            ('Event','$t_{\\nu}$-$t_{trigger}$ (s)','RA','Dec',
                            '\sigma (90\%)','E_{reco} (GeV)',
                            'p-value','In 90\% Contour')]
                    else:
                        event_table+=[
                            ('$t_{\\nu}$-$t_{trigger}$ (s)','RA','Dec',
                            '\sigma (90\%)','E_{reco} (GeV)',
                            'p-value','In 90\% Contour')]
                else:
                    with_p = False
                    event_table+=[
                        ('$t_{\\nu}$-$t_{trigger}$ (s)','RA','Dec',
                        '\sigma (90\%)','E_{reco} (GeV)',
                        'In 90\% Contour')]
                i = 0
                for event in self.analysis.coincident_events:
                    if with_p:
                        if len(self.analysis.coincident_events)>100:
                            if event['in_contour']:
                                event_table+=[
                                ('{:.0f}'.format((event['time']-self.source['trigger_mjd'])*86400.),
                                "{:3.2f}\degree".format(np.rad2deg(event['ra'])),
                                '{:3.2f}\degree'.format(np.rad2deg(event['dec'])),
                                "{:3.2f}\degree".format(np.rad2deg(event["sigma"]*self.analysis._angScale)),
                                "{:.2E}".format(10**event['logE']),
                                '{:.4f}'.format(event['pvalue']),
                                str(bool(event['in_contour']))
                                )]
                        else: 
                            event_table+=[
                            ('{}'.format(i+1),
                            '{:.0f}'.format((event['time']-self.source['trigger_mjd'])*86400.),
                            "{:3.2f}\degree".format(np.rad2deg(event['ra'])),
                            '{:3.2f}\degree'.format(np.rad2deg(event['dec'])),
                            "{:3.2f}\degree".format(np.rad2deg(event["sigma"]*self.analysis._angScale)),
                            "{:.2E}".format(10**event['logE']),
                            '{:.4f}'.format(event['pvalue']),
                            str(bool(event['in_contour']))
                            )]
                            i+=1
                    else:
                        if len(self.analysis.coincident_events)>100:
                            if event['in_contour']:
                                event_table+=[
                                ('{:.0f}'.format((event['time']-self.source['trigger_mjd'])*86400.),
                                "{:3.2f}\degree".format(np.rad2deg(event['ra'])),
                                '{:3.2f}\degree'.format(np.rad2deg(event['dec'])),
                                "{:3.2f}\degree".format(np.rad2deg(event["sigma"]*self.analysis._angScale)),
                                "{:.2E}".format(10**event['logE']),
                                str(bool(event['in_contour']))
                                )]
                        else: 
                            event_table+=[
                            ('{:.0f}'.format((event['time']-self.source['trigger_mjd'])*86400.),
                            "{:3.2f}\degree".format(np.rad2deg(event['ra'])),
                            '{:3.2f}\degree'.format(np.rad2deg(event['dec'])),
                            "{:3.2f}\degree".format(np.rad2deg(event["sigma"]*self.analysis._angScale)),
                            "{:.2E}".format(10**event['logE']),
                            str(bool(event['in_contour']))
                            )]
            self.write_table(f,"event",[],event_table)
        else:
            f.write(r"\newcommand{\event}{[None]}")

    def make_new_command(self, f, name, command):
        """Define a new command, in LaTeX, for the report
        
        Parameters
        -----------
        f: I/O file object
            Opened file that is being written to
        name: str
            Name of the command being defined
        command: str
            Value of the command being defined
        """
        f.write(r"\newcommand{"+"\\"+name+"}{"+
                command+"}\n")

    def generate_report(self):
        """Generate the report! This function does most of the heavy lifting.
        Starts making the source LaTeX file (named r.tex), queries runs,
        makes detector status plots, includes all plots for the analysis,
        makes all tables for analysis parameters, runs, events, etc.
        """

        report_fname = os.path.join(self.dirname, "r.tex")

        self.run_table = self.query_db_runs(self.time_window)
        now = datetime.datetime.utcnow()
        self.ontime = {}
        self.ontime['type'] = 'database'
        self.ontime['stream'] = self.stream
        self.ontime['runkey'] = 'run_id'
        self.ontime['time_start'] = Time(
            self.run_table[0]['start'],
            format='iso',
            scale='utc',
            precision=0).iso
        self.ontime['time_stop'] =Time(
            self.run_table[-1]['stop'],
            format='iso',
            scale='utc',
            precision=0).iso
        self.ontime['time_query']=now.strftime('%Y-%m-%d %H:%M:%S')

        self.query_events=icecube.realtime_tools.live.get_events(
            self.ontime['stream'],
            self.ontime['time_start'],
            self.ontime['time_stop'])

        self.widewindow = self.ontime_table(self.query_events)
        try:
            self.widewindow['t']=Time(
                list(self.widewindow['eventtime']),
                scale='utc',
                format='iso')
            for run in self.run_table:
                run['gfu_counts'] = (self.widewindow['run_id']==run['run_number']).sum()
        except:
            print("Old years of data have different database keys")
            for run in self.run_table:
                run['gfu_counts'] = 0. 

        s = self.source

        dataset = Datasets[self.analysis._dataset]
        # Make rate plots and put them in analysis directory
        make_rate_plots(
            self.time_window,
            self.run_table,
            self.query_events,
            self.dirname,
            season=self.stream)

        with open(report_fname, 'wt') as f:

            d1 = s["start_iso"].split()[0]
            d2 = s["stop_iso"].split()[0]

            if d1 == d2:
                obsdate = d1
            else:
                obsdate = d1+' -- '+d2

            for name, command in [('sourcename', self.analysis.name),
                                  ('analysisid', self.analysisid),
                                  ('reportdate', datetime.date.today().strftime("%Y-%m-%d")),
                                  ('obsdate', obsdate)]:
                self.make_new_command(f, name, command)

            for plot_name, plot_path in self._figure_list:
                f.write(
                    r"\newcommand{"+"\\"+plot_name+"}{"+ "\\includegraphics[width=0.8\\textwidth]" +
                    "{" + self.dirname+"/"+
                    plot_path + "}" +
                    "}\n"
                )

            if self.source["stop_mjd"] > 60276.86875:
                muonfilter_str = "MuonFilter_23_plot.png"
                l2filter_str="OnlineL2Filter_23_plot.png"
            else: 
                muonfilter_str="MuonFilter_13_plot.png"
                l2filter_str="OnlineL2Filter_17_plot.png"

            for plot_name, plot_path in [("gfurate", "GFU_rate_plot.png"),
                                         ('skymap', self.analysisid + "unblinded_skymap.png"),
                                         ('skymapzoom', self.analysisid + "unblinded_skymap_zoom.png"),
                                         ('limitdNdE', "central_90_dNdE.png"),
                                         ('muonfilter', muonfilter_str),
                                         ("Lfilter", l2filter_str),
                                         ("badnessplot", "badness_plot.png"),
                                         ("multiplicity", "IN_ICE_SIMPLE_MULTIPLICITY_plot.png")]:
                if os.path.isfile(self.dirname + '/' + plot_path):
                    f.write(
                        r"\newcommand{"+"\\"+plot_name+"}{"+ self.dirname+"/"+
                        plot_path +
                        "}\n"
                    )

            for plot_name, plot_path, else_message in [('tsd', 'TS_distribution.png', 'No background TS distribution'),
                                                       ('upperlim', 'upper_limit_distribution.png', "No upper limit calculation")]:
                if os.path.isfile(self.dirname + '/' + plot_path):
                    f.write(
                        r"\newcommand{"+"\\"+plot_name+"}{"+ "\\includegraphics[width=0.8\\textwidth]" +
                        "{" + self.dirname+"/"+
                        plot_path + "}" +
                        "}\n"
                    )
                else:
                    f.write(
                        r"\newcommand{"+"\\"+plot_name+"}{"+ 
                            else_message + 
                        "}\n"
                    )

            if self.source_type == "PS":
                location_string = f'{np.rad2deg(self.analysis.ra):.2f}, {np.rad2deg(self.analysis.dec):+.2f}'
                location_entry = "Source location (deg.)"
            else:
                location_string = f'{self.analysis.skymap_path}'
                location_entry = "Skymap path"
            
            self.write_table(
                f,"sourcetable", [],[
                    ("Source Name", self.analysis.name),
                    (location_entry, location_string.replace('_', '\_')),
                    ("Trigger Time", "{} (MJD={:12.6f})".format(
                        s["trigger_iso"], s["trigger_mjd"])),
                    ("Start Time", "{} (Trigger{:+1.1f}s)".format(
                        s["start_iso"], (self.time_window[0].mjd-s['trigger_mjd'])*86400.)),
                    ("Stop Time", "{} (Trigger{:+1.1f}s)".format(
                        s["stop_iso"], (self.time_window[1].mjd-s['trigger_mjd'])*86400.)),
                    ("Time Window",r"{:1.1f}s".format(s["realtime"])),
                ]
            )

            self.write_table(f,"skylabtable",[],[
                ("Skylab Version", skylab.__version__),
                ("IceTray Path", str(icetray.__path__).replace('_', '\_')),
                ("Created by", expanduser('~')[6:]),
                ("Dataset Used", str(dataset.subdir).replace('_',' ')),
                ("Dataset details", str(dataset.name)[:80]),
                ("", str(dataset.name)[80:])
            ])

            r1=[]
            r2=[]
            livetime = 0

            for run in self.run_table:
                r1.append((
                    run["run_number"],
                    run["start"].split(".")[0],
                    run["stop"].split(".")[0],
                    run["duration"],"{:1.1f}s".format(run["livetime"])
                ))

                livetime += run['livetime']
            
            self.write_table(
                f,
                "runtimetable",
                ["Run","Start Time","Stop Time","Duration","Livetime"],
                r1
            )

            if 'status' in self.run_table[0]:
                for run in self.run_table:
                    r2.append((run['run_number'],run['status'],run['lightmode'],run['filter_mode'],
                               run['run_mode'],run["OK"],run["gfu_counts"]))
                self.write_table(
                    f,
                    "runstatustable",
                    ["Run","Status","Light","Filter Mode","Run Mode","OK","GFU" ],
                    r2
                )
            else:
                self.write_table(f,"runstatustable",[],[])

            f.write(r"\newcommand{\livetime}{"+'{:0,.1f}'.format(livetime)+"}\n")

            if self.ontime['type']=='database':
                self.write_table(
                    f, "ontimetable", [], [
                        ("Access Method",self.ontime['type']),
                        ("Stream", r"\texttt{"+self.ontime['stream']+"}"),
                        ("Query Time",self.ontime['time_query']),
                        ("Start Time",self.ontime['time_start']),
                        ("Stop Time", self.ontime['time_stop']),
                    ]
                )

            self.make_coinc_events_table(f)
            
            if self.analysis._float_index:
                if self.analysis.ts > 0.:
                    self.write_table(
                        f,
                        "results",
                        [],
                        [("$n_s$", "{:1.3f}".format(self.analysis.ns)),
                        ("$TS$", "{:1.3f}".format(self.analysis.ts)),
                        ("$\gamma$", f"{self.analysis.gamma:.2f}"),
                        ("$p-value$", "{:1.4f}".format(self.analysis.p)),
                        ("best-fit RA", "{:3.2f}\degree".format(np.rad2deg(self.analysis.skymap_fit_ra))),
                        ("best-fit dec", "{:3.2f}\degree".format(np.rad2deg(self.analysis.skymap_fit_dec)))]
                    )
                else: 
                    self.write_table(
                        f,
                        "results",
                        [],
                        [("$n_s$", "{:1.3f}".format(self.analysis.ns)),
                        ("$TS$", "{:1.3f}".format(self.analysis.ts)),
                        ("$\gamma$", f"{self.analysis.gamma:.2f}"),
                        ("$p-value$", "{:1.4f}".format(self.analysis.p))]
                    )
            else:
                self.write_table(
                    f,
                    "results",
                    [],
                    [("$n_s$","{:1.3f}".format(self.analysis.ns)),
                     ("$TS$","{:1.3f}".format(self.analysis.ts)),
                     ("$p-value$","{:1.4f}".format(self.analysis.p))]
                )

    def make_pdf(self):
        """Makes the PDF compiled version of the LaTeX report
        using pdflatex, and saves to the analysis folder
        """
        # symlink main report tex file
        reportfname = self.analysisid + "_report.tex"
        reportpath = os.path.join(self.dirname, reportfname)
        reportsrc = self.get_report_source()
        
        if os.path.exists(reportpath):
            os.unlink(reportpath)

        os.symlink(reportsrc, reportpath)

        # get environment variables
        env = dict(os.environ)
        subprocess.call(
            ['pdflatex', '-interaction=batchmode',
             '-output-directory=%s' % self.dirname,
             self.dirname + '/' + self.analysisid + "_report.tex"],
            #cwd=self.dirname,
            env = env,
        )

class FastResponseReport(ReportGenerator):
    """Generates a general FastResponse Report for a point source.
    Includes additional plots of the 1D LLH ns scan
    """
    _figure_list = [('nsscan', 'llh_ns_scan.png')]

    def __init__(self, analysis):
        super().__init__(analysis)

    def generate_report(self):
        super().generate_report()

class GravitationalWaveReport(ReportGenerator):
    """Generates a FastResponse Report for a GW.
    Includes additional plot of skymap PDF versus dec, 
    with PS sensitivity overlaid
    """
    _figure_list = [('decPDF', 'decPDF.png')]

    def get_report_source(self):
        """Find base path of the general report TeX document.
        Looks for the directory location fast_response, 
        returns fast_response/latex/report_gw.tex
        """
        base = os.path.dirname(fast_response.__file__)
        return os.path.join(base, 'latex/report_gw.tex')

class AlertReport(ReportGenerator):
    """Generates a FastResponse Report for an Alert event."""
    _figure_list = []

    def get_report_source(self):
        """Find base path of the general report TeX document.
        Looks for the directory location fast_response, 
        returns fast_response/latex/report_alert.tex
        """
        base = os.path.dirname(fast_response.__file__)
        return os.path.join(base, 'latex/report_alert.tex')
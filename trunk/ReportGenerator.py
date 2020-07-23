'''Class for generating reports for Fast Reponse and GW
   followup analyses. A lot of code has been taken from Fast Reponse
   Analysis. 

   Author: Kevin Meagher, Raamis Hussain & Alex Pizzuto
   Date:   Mar 27, 2019
   Modified: April, 2019
'''

import numpy as np
import pandas as pd
import os,json,hashlib,warnings,datetime,types,marshal
import re
import logging as log
import urllib.request, urllib.parse, urllib.error,urllib.request,urllib.error,urllib.parse
import icecube.realtime_tools.live
import subprocess
import sys
import skylab
from skylab.datasets import Datasets
from os.path import expanduser

from icecube             import astro
from icecube             import icetray
from astropy.time        import Time,TimeDelta
from astropy.coordinates import Angle
from astropy             import units
from make_ontime_plots   import make_rate_plots

log.basicConfig(level=log.ERROR)
mpl_logger = log.getLogger('matplotlib') 
mpl_logger.setLevel(log.ERROR) 
 
def interptime(tstr):
    if tstr is None:
        t= None
    elif tstr.endswith('sec'):
        t = TimeDelta(float(tstr[:-3]),format='sec')
    elif tstr.endswith('s'):
        t = TimeDelta(float(tstr[:-1]),format='sec')
    elif tstr.endswith('min'):
        t = TimeDelta(float(tstr[:-3])*60,format='sec')
    elif tstr.endswith('hr'):
        t = TimeDelta(float(tstr[:-2])*3600,format='sec')
    elif tstr.startswith('mjd'):
        t = Time(float(tstr[3:]),format='mjd')
    elif tstr.endswith('d'):
        t = TimeDelta(float(tstr[:-1]),format='jd')
    else:
        t = Time(tstr,format='iso')
    return t
 
def get_times(trig,sta,sto):

    trigger = interptime(trig)

    if type(trigger)!=Time:
        trigger=None

    start = interptime(sta)

    if type(trigger)==Time and type(start)==Time:
        pass
    elif type(trigger)==Time and start is None:
        start = trigger
    elif type(start)==Time and trigger is None:
        trigger = start
    elif type(start)==TimeDelta and trigger is not None:
        start = trigger+start
    else:
        raise Exception("either --trigger or --start must be a absolute time")

    stop=interptime(sto)
    if type(stop)!=Time:
        stop = trigger+stop

    return trigger,start,stop

def valid_location(ra, dec):
    r'''For a point source, check if location
    is valid
    Parameters:
    -----------
    ra: float
        Right ascension in radians
    dec: float
        Declination in radians
    Returns:
    --------
    True or False: Bool
        Whether or not the location is valid
    '''
    if (ra is None) or (dec is None):
        return False
    else:
        if (ra < 0.) or (ra > 2*np.pi):
            return False
        elif (dec < -np.pi / 2.) or (dec > np.pi/2.):
            return False
        else:
            return True


class ReportGenerator(object):

    def __init__(self,name,trigger,start,stop,ts,ns,source_type,id,**kwargs):
        r""" Constructor

        Parameters
        ----------
        
        """ 
        self.source = {}
        self.source['name'] = name
        self.source['trigger'] = trigger
        self.source['start'] = start
        self.source['stop'] = stop
        self.source['ts'] = ts
        self.source['ns'] = ns
        self.source['source_type'] = source_type

        self.source['pvalue'] = kwargs.pop('p', None)
        self.source['sigma'] = kwargs.pop('sigma', None)

        if self.source['source_type'] == "PS":
            self.source['ra'] = kwargs.pop("ra", None)
            self.source['dec'] = kwargs.pop("dec", None)
            self.source['extension'] = kwargs.pop("extension", None)
            self.skymap = None
            if not valid_location(self.source['ra'], self.source['dec']):
                print("Invalid location for a point source")
                sys.exit()
            else:
                self.source["ra_deg"]=np.rad2deg(self.source['ra'])
                self.source["dec_deg"]=np.rad2deg(self.source['dec'])
            if self.source['extension'] is not None:
                self.source['extension'] = np.rad2deg(self.source['extension'])
        elif self.source['source_type'] == "Prior":
            self.source['skymap'] = kwargs.pop("skymap", None)
            self.source['ra'], self.source['dec'] = None, None
            if self.source['skymap'] is None:
                print("Spatial Prior Analyses require skymaps")
                sys.exit()
        else:
            print("Invalid source type")
            sys.exit()

        try:
            assert stop > start
        except:
            print("Stop Time must come after Start time. Exiting")
            sys.exit()

        self.source["time_trigger_iso"] = Time(self.source['trigger'], format='mjd').iso
        self.source["time_start_iso"] = Time(self.source['start'], format='mjd').iso
        self.source["time_stop_iso"] = Time(self.source['stop'], format='mjd').iso
        self.source["time_trigger_mjd"] = self.source['trigger']
        self.source["time_start_mjd"] = self.source['start']
        self.source["time_stop_mjd"] = self.source['stop']
        self.source["realtime"]= (stop - start) * 86400. #Convert days to seconds
        if self.source["time_stop_mjd"] > 58309.74:
            self.stream = 'neutrino'
        elif self.source["time_stop_mjd"] > 57891.17:
            self.stream = 'neutrino17'
        else:
            self.stream = 'neutrino16'

        self.time_window=(Time(self.source['start'], format='mjd'), Time(self.source['stop'], format='mjd'))
        self.trigger = Time(self.source['trigger'], format='mjd')
        self.analysisid = id
        self.dirname = kwargs.pop('analysispath','./')
        self.events = kwargs.pop('coincident_events',None)

    def write_table(self,file,name,header,table,prefix=""):
        """
        write a table in latex format for the generated reports
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

    def query_db_runs(self,time_window):
        # first query the rundata base with the times of the analysis
        # plus 8 hours on either side
        run_url = 'https://live.icecube.wisc.edu/run_info/'
        run_query = {'user':'icecube', 'pass':'skua',
                     'start':(Time(time_window[0],precision=0)
                              - TimeDelta(8*3600,format='sec')).iso,
                     'stop': (Time(time_window[1],precision=0)
                              + TimeDelta(8*3600,format='sec')).iso}
        now = Time(datetime.datetime.now()+datetime.timedelta(hours=6),scale='utc',precision=0)
        req_data = urllib.parse.urlencode(run_query).encode("utf-8")
        req = urllib.request.Request(run_url)
        with urllib.request.urlopen(req, data=req_data, timeout=500) as fi:
            run_table = json.load(fi)
        #run_table= json.loads(
        #        urllib.request.urlopen(urllib.request.Request(
        #            run_url,urllib.parse.urlencode(run_query)),timeout=500).read())

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
            if run['livetime'] > 86400. * 2:
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
        newdict=[]
        for event in query_dict:
            newevent = event['value']['data']
            for key,val in list(newevent['reco']['splinempe'].items()): 
                newevent['splinempe_'+key]=val
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

    def generate_gw_report(self):

        report_fname = os.path.join(self.dirname,"r.tex")

        self.run_table = self.query_db_runs(self.time_window)
        now = datetime.datetime.utcnow()
        self.ontime = {}
        self.ontime['type']='database'
        self.ontime['stream'] = self.stream
        self.ontime['runkey']='run_id'
        self.ontime['time_start']=Time(self.run_table[0]['start'],format='iso',scale='utc',precision=0).iso
        self.ontime['time_stop'] =Time(self.run_table[-1]['stop'],format='iso',scale='utc',precision=0).iso
        self.ontime['time_query']=now.strftime('%Y-%m-%d %H:%M:%S')

        self.query_events=icecube.realtime_tools.live.get_events(
            self.ontime['stream'], self.ontime['time_start'],self.ontime['time_stop'])

        self.widewindow = self.ontime_table(self.query_events)
        self.widewindow['t']=Time(list(self.widewindow['eventtime']),scale='utc',format='iso')

        for run in self.run_table:
            run['gfu_counts'] = (self.widewindow['run_id']==run['run_number']).sum()

        s=self.source

        # Make rate plots and put them in analysis directory
        make_rate_plots(self.time_window,self.run_table,self.query_events,self.dirname)

        with open(report_fname,'wt') as f:

            d1=s["time_start_iso"].split()[0]
            d2=s["time_stop_iso"].split()[0]

            if d1==d2:
                obsdate = d1
            else:
                obsdate = d1+' -- '+d2

            f.write(
                r"\newcommand{"+"\\"+"sourcename"+"}{"+
                s["name"]+"}\n")

            f.write(r"\newcommand{\analysisid}{"+self.analysisid+"}\n")
            
            f.write(
                r"\newcommand{"+"\\"+"gfurate"+"}{"+ self.dirname+"/"+
                "GFU_rate_plot.png"+
                "}\n")

            f.write(
                r"\newcommand{"+"\\"+"reportdate"+"}{"+
                datetime.date.today().strftime("%Y-%m-%d")+
                "}\n")

            f.write(
                r"\newcommand{"+"\\"+"skymap"+"}{"+ self.dirname+"/"+
                self.analysisid + "unblinded_skymap.png" +
                "}\n")

            f.write(
                r"\newcommand{"+"\\"+"muonfilter"+"}{"+ self.dirname+"/"+
                "MuonFilter_13_plot.png"+
                "}\n")
       
            f.write(
                r"\newcommand{"+"\\"+"Lfilter"+"}{"+ self.dirname+"/"+
                "OnlineL2Filter_17_plot.png"+
                "}\n")
            
            f.write(
                r"\newcommand{"+"\\"+"badness"+"}{"+ self.dirname+"/"+
                "badness_plot.png"+
                "}\n")

            f.write(
                r"\newcommand{"+"\\"+"obsdate"+"}{"+
                obsdate+"}\n")

            f.write(
                r"\newcommand{"+"\\"+"multiplicity"+"}{"+ self.dirname+"/"+
                "IN_ICE_SIMPLE_MULTIPLICITY_plot.png"+
                "}\n")

            sf_fname = os.path.join(self.dirname,self.analysisid+"_SurvivalFunction.pdf")
            if os.path.exists(sf_fname):
                f.write(r"\newcommand{"+"\\"+"survivialfunctionplot"+"}{"+
                        r"\includegraphics[width=\textwidth]{\analysisid_SurvivalFunction}"+
                        "}\n")
            else:
                f.write(r"\newcommand{"+"\\"+"survivialfunctionplot"+"}{}\n")

            dist_fname = os.path.join(self.dirname,self.analysisid+"_BackgroundPDF.pdf")
            if os.path.exists(sf_fname):
                f.write(r"\newcommand{"+"\\"+"backgroundpdfplot"+"}{"+
                        r"\includegraphics[width=\textwidth]{\analysisid_BackgroundPDF}"+
                        "}\n")
            else:
                f.write(r"\newcommand{"+"\\"+"backgroundpdfplot"+"}{}\n")

            self.write_table(f,"sourcetable",[],[
                ("Source Name",s["name"]),
                ("Trigger Time","{} (MJD={:12.6f})".format(s["time_trigger_iso"],s["time_trigger_mjd"])),
                ("Start Time", "{} (Trigger{:+1.1f}s)".format(s["time_start_iso"],
                                                              (self.time_window[0]-self.trigger).sec)),
                ("Stop Time", "{} (Trigger{:+1.1f}s)".format(s["time_stop_iso"],
                                                             (self.time_window[1]-self.trigger).sec)),
                ("Time Window",r"{:1.1f}s".format(s["realtime"])),
            ])

            r1=[]
            r2=[]
            livetime = 0

            for run in self.run_table:
                r1.append((run["run_number"],run["start"].split(".")[0],run["stop"].split(".")[0],run["duration"],"{:1.1f}s".format(run["livetime"])))

                livetime += run['livetime']
            self.write_table(f,"runtimetable",["Run","Start Time","Stop Time","Duration","Livetime"],r1,)

            if 'status' in self.run_table[0]:
                for run in self.run_table:
                    r2.append((run['run_number'],run['status'],run['lightmode'],run['filter_mode'],
                               run['run_mode'],run["OK"],run["gfu_counts"]))
                self.write_table(f,"runstatustable",["Run","Status","Light","Filter Mode","Run Mode","OK","GFU" ],r2,)
            else:
                self.write_table(f,"runstatustable",[],[])

            f.write(r"\newcommand{\livetime}{"+'{:0,.1f}'.format(livetime)+"}\n")

            if self.ontime['type']=='database':
                self.write_table(f,"ontimetable",[],[("Access Method",self.ontime['type']),
                                                ("Stream", r"\texttt{"+self.ontime['stream']+"}"),
                                                ("Query Time",self.ontime['time_query']),
                                                ("Start Time",self.ontime['time_start']),
                                                ("Stop Time", self.ontime['time_stop']),
                                                ])


            event_table = []
            if self.events is not None:
                for event in self.events:
                    event_table+=[
                        ("Run",'{}'.format(event['run'])),
                        ("Event",'{}'.format(event['event'])),
                        ("Time","{}".format(
                            event['time'])),
                        ("Right Ascension","{:3.2f}\degree"
                         .format(np.rad2deg(event['ra']))),
                        ("Declination","{:3.2f}\degree"
                         .format(np.rad2deg(event['dec']))),
                        ("Angular Uncertainty (90\%)","{:3.2f}\degree"
                            .format(np.rad2deg(event["sigma"]*2.145966))),
                        ("Reconstructed Energy (GeV)","{:2.2f}"
                            .format(10**event['logE'])),
                        None,
                    ]

                self.write_table(f,"event",[],event_table)
            else:
                f.write(r"\newcommand{\event}{[None]}")

            llh_table = []
            for i in range(s['ts'].size):
                llh_table+=[
                    ("$TS$",'{:1.3f}'.format(s['ts'][i])),
                    ("$n_s$",'{:1.3f}'.format(s['ns'][i])),
                    ("P value","{:1.4f}".format(s['pvalue'][i])),
                    ("$\gamma$","{:1.3f}".format(s['gamma'][i]))
                ]

            self.write_table(f,"results",[],llh_table)
            # self.write_table(f,"results",[],[
            #     ("$n_s$","{:1.3f}".format(s['ns'])),
            #     ("$TS$","{:1.3f}".format(s['ts'])),
            #     ("$P value$","{:1.4f}".format(s['pvalue'])),
            #     ("$\gamma$","{:1.3f}".format(s['gamma']))
            # ])           


        # symlink main report tex file
        reportfname = self.analysisid+"_report.tex"
        reportpath = os.path.join(self.dirname,reportfname)
        reportsrc = os.path.join(os.environ["I3_BUILD"],'fast_response','resources','latex','report_gw.tex')
        if os.path.exists(reportpath):
            os.unlink(reportpath)

        os.symlink(reportsrc,reportpath)

        # symlink directory of needed sty files
        styledname = 'sty'
        styledpath = os.path.join(self.dirname, styledname)
        styledsrc  = os.path.join(os.environ["I3_BUILD"],'fast_response','resources','latex','sty')
        if os.path.exists(styledpath):
            os.unlink(styledpath)
        # but only if source exists (user decision, only needed on some systems)
        if os.path.exists(styledsrc):
            os.symlink(styledsrc,styledpath)

    def generate_report(self):

        report_fname = os.path.join(self.dirname,"r.tex")

        self.run_table = self.query_db_runs(self.time_window)
        now = datetime.datetime.utcnow()
        self.ontime = {}
        self.ontime['type']='database'
        self.ontime['stream'] = self.stream
        self.ontime['runkey']='run_id'
        self.ontime['time_start']=Time(self.run_table[0]['start'],format='iso',scale='utc',precision=0).iso
        self.ontime['time_stop'] =Time(self.run_table[-1]['stop'],format='iso',scale='utc',precision=0).iso
        self.ontime['time_query']=now.strftime('%Y-%m-%d %H:%M:%S')

        self.query_events=icecube.realtime_tools.live.get_events(
            self.ontime['stream'], self.ontime['time_start'],self.ontime['time_stop'])

        self.widewindow = self.ontime_table(self.query_events)
        try:
            self.widewindow['t']=Time(list(self.widewindow['eventtime']),scale='utc',format='iso')
            for run in self.run_table:
                run['gfu_counts'] = (self.widewindow['run_id']==run['run_number']).sum()
        except:
            print("Old years of data have different database keys")
            for run in self.run_table:
                run['gfu_counts'] = 0. 

        s=self.source

        dataset = Datasets['GFUOnline']
        # Make rate plots and put them in analysis directory
        make_rate_plots(self.time_window,self.run_table,self.query_events,self.dirname, season=self.stream)

        with open(report_fname,'wt') as f:

            d1=s["time_start_iso"].split()[0]
            d2=s["time_stop_iso"].split()[0]

            if d1==d2:
                obsdate = d1
            else:
                obsdate = d1+' -- '+d2

            f.write(
                r"\newcommand{"+"\\"+"sourcename"+"}{"+
                s["name"]+"}\n")

            f.write(r"\newcommand{\analysisid}{"+self.analysisid+"}\n")
            
            f.write(
                r"\newcommand{"+"\\"+"gfurate"+"}{"+ self.dirname+"/"+
                "GFU_rate_plot.png"+
                "}\n")

            f.write(
                r"\newcommand{"+"\\"+"reportdate"+"}{"+
                datetime.date.today().strftime("%Y-%m-%d")+
                "}\n")

            f.write(
                r"\newcommand{"+"\\"+"skymap"+"}{"+ self.dirname+"/"+
                self.analysisid + "unblinded_skymap.png" +
                "}\n")
            
            f.write(
                r"\newcommand{"+"\\"+"skymapzoom"+"}{"+ self.dirname+"/"+
                self.analysisid + "unblinded_skymap_zoom.png" +
                "}\n")

            f.write(
                r"\newcommand{"+"\\"+"limitdNdE"+"}{"+ self.dirname+"/"+
                "central_90_dNdE.png" +
                "}\n"
            )

            f.write(
                r"\newcommand{"+"\\"+"muonfilter"+"}{"+ self.dirname+"/"+
                "MuonFilter_13_plot.png"+
                "}\n")
       
            f.write(
                r"\newcommand{"+"\\"+"Lfilter"+"}{"+ self.dirname+"/"+
                "OnlineL2Filter_17_plot.png"+
                "}\n")
            
            f.write(
                r"\newcommand{"+"\\"+"badnessplot"+"}{"+ self.dirname+"/"+
                "badness_plot.png"+
                "}\n")
            
            if os.path.isfile(self.dirname + '/TS_distribution.png'):
                f.write(
                    r"\newcommand{"+"\\"+"tsd"+"}{"+ "\\includegraphics[width=0.9\\textwidth]" +
                    "{" + self.dirname+"/"+
                    "TS_distribution.png"+ "}" +
                    "}\n")
            else:
                f.write(
                    r"\newcommand{"+"\\"+"tsd"+"}{"+ 
                        "No background TS distribution" + 
                    "}\n"
                )

            if os.path.isfile(self.dirname + '/upper_limit_distribution.png'):
                f.write(
                r"\newcommand{"+"\\"+"upperlim"+"}{"+ "\\includegraphics[width=0.9\\textwidth]" +
                    "{" + self.dirname+"/"+
                    "upper_limit_distribution.png"+ "}" +
                    "}\n"
                )
            else:
                f.write(
                    r"\newcommand{"+"\\"+"upperlim"+"}{"+
                        "No upper limit calculation" +
                    "}\n"
                )

            if os.path.isfile(self.dirname + '/llh_ns_scan.png'):
                f.write(
                r"\newcommand{"+"\\"+"nsscan"+"}{"+ "\\includegraphics[width=0.9\\textwidth]" +
                    "{" + self.dirname+"/"+
                    "llh_ns_scan.png"+ "}" +
                    "}\n" 
                )
            else:
                f.write(
                    r"\newcommand{"+"\\"+"nsscan"+"}{"+
                        "No llh vs. ns scan" +
                    "}\n"
                )

            f.write(
                r"\newcommand{"+"\\"+"obsdate"+"}{"+
                obsdate+"}\n")

            f.write(
                r"\newcommand{"+"\\"+"multiplicity"+"}{"+ self.dirname+"/"+
                "IN_ICE_SIMPLE_MULTIPLICITY_plot.png"+
                "}\n")

            sf_fname = os.path.join(self.dirname,self.analysisid+"_SurvivalFunction.pdf")
            if os.path.exists(sf_fname):
                f.write(r"\newcommand{"+"\\"+"survivialfunctionplot"+"}{"+
                        r"\includegraphics[width=\textwidth]{\analysisid_SurvivalFunction}"+
                        "}\n")
            else:
                f.write(r"\newcommand{"+"\\"+"survivialfunctionplot"+"}{}\n")

            dist_fname = os.path.join(self.dirname,self.analysisid+"_BackgroundPDF.pdf")
            if os.path.exists(sf_fname):
                f.write(r"\newcommand{"+"\\"+"backgroundpdfplot"+"}{"+
                        r"\includegraphics[width=\textwidth]{\analysisid_BackgroundPDF}"+
                        "}\n")
            else:
                f.write(r"\newcommand{"+"\\"+"backgroundpdfplot"+"}{}\n")

            self.write_table(f,"sourcetable",[],[
                ("Source Name",s["name"]),
                ("Trigger Time","{} (MJD={:12.6f})".format(s["time_trigger_iso"],s["time_trigger_mjd"])),
                ("Start Time", "{} (Trigger{:+1.1f}s)".format(s["time_start_iso"],
                                                              (self.time_window[0]-self.trigger).sec)),
                ("Stop Time", "{} (Trigger{:+1.1f}s)".format(s["time_stop_iso"],
                                                             (self.time_window[1]-self.trigger).sec)),
                ("Time Window",r"{:1.1f}s".format(s["realtime"])),
            ])

            self.write_table(f,"skylabtable",[],[
                ("Skylab Version",skylab.__version__),
                ("IceTray Path",str(icetray.__path__).replace('_', '\_')),
                ("Created by", expanduser('~')[6:]),
                ("Dataset Used",str(dataset.subdir).replace('_',' ')),
                ("Dataset details", str(dataset.name)[:80]),
                ("", str(dataset.name)[80:])
            ])

            r1=[]
            r2=[]
            livetime = 0

            for run in self.run_table:
                r1.append((run["run_number"],run["start"].split(".")[0],run["stop"].split(".")[0],run["duration"],"{:1.1f}s".format(run["livetime"])))

                livetime += run['livetime']
            self.write_table(f,"runtimetable",["Run","Start Time","Stop Time","Duration","Livetime"],r1,)

            if 'status' in self.run_table[0]:
                for run in self.run_table:
                    r2.append((run['run_number'],run['status'],run['lightmode'],run['filter_mode'],
                               run['run_mode'],run["OK"],run["gfu_counts"]))
                self.write_table(f,"runstatustable",["Run","Status","Light","Filter Mode","Run Mode","OK","GFU" ],r2,)
            else:
                self.write_table(f,"runstatustable",[],[])

            f.write(r"\newcommand{\livetime}{"+'{:0,.1f}'.format(livetime)+"}\n")

            if self.ontime['type']=='database':
                self.write_table(f,"ontimetable",[],[("Access Method",self.ontime['type']),
                                                ("Stream", r"\texttt{"+self.ontime['stream']+"}"),
                                                ("Query Time",self.ontime['time_query']),
                                                ("Start Time",self.ontime['time_start']),
                                                ("Stop Time", self.ontime['time_stop']),
                                                ])


            event_table = []
            if self.events is not None:
                for event in self.events:
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

                self.write_table(f,"event",[],event_table)
            else:
                f.write(r"\newcommand{\event}{[None]}")

            #llh_table = []
            #for i in range(s['ts'].size):
            #    llh_table+=[
            #        ("$TS$",'{:1.3f}'.format(s['ts'][i])),
            #        ("$n_s$",'{:1.3f}'.format(s['ns'][i])),
            #        ("P value","{:1.4f}".format(s['pvalue'][i])),
            #        ("$\gamma$","{:1.3f}".format(s['gamma'][i]))
            #    ]

            #self.write_table(f,"results",[],llh_table)
            self.write_table(f,"results",[],[
                 ("$n_s$","{:1.3f}".format(s['ns'])),
                 ("$TS$","{:1.3f}".format(s['ts'])),
                 ("$p-value$","{:1.4f}".format(s['pvalue']))
            ])           


        # symlink main report tex file
        reportfname = self.analysisid+"_report.tex"
        reportpath = os.path.join(self.dirname,reportfname)
        #reportsrc = os.path.join(os.environ["I3_BUILD"],'fast_response','resources','latex','report_skylab.tex')
        try:
            reportsrc = os.path.join(os.environ["FAST_RESPONSE_SCRIPTS"], 'latex', 'report_skylab.tex')
        except KeyError:
            reportsrc = os.path.join(os.getcwd(), 'latex', 'report_skylab.tex')
        if os.path.exists(reportpath):
            os.unlink(reportpath)

        os.symlink(reportsrc,reportpath)

        # symlink directory of needed sty files
        styledname = 'sty'
        styledpath = os.path.join(self.dirname, styledname)
        #styledsrc  = os.path.join(os.environ["I3_BUILD"],'fast_response','resources','latex','sty')
        try:
            styledsrc = os.path.join(os.environ["FAST_RESPONSE_SCRIPTS"], 'latex', 'sty')
        except KeyError:
            styledsrc = os.path.join(os.getcwd(), 'latex', 'sty')
        if os.path.exists(styledpath):
            os.unlink(styledpath)
        # but only if source exists (user decision, only needed on some systems)
        if os.path.exists(styledsrc):
            os.symlink(styledsrc,styledpath)

    def make_pdf(self):
        # get environment variables
        env = dict(os.environ)
        # add location for sty files
        env['TEXINPUTS'] = './:./sty/:'
        subprocess.call(['pdflatex','-interaction=batchmode','-output-directory=%s' % self.dirname, 
                        self.dirname+'/'+self.analysisid+"_report.tex"],
                        #cwd=self.dirname,
                        env = env,
                       )

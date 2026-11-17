'''Class for generating reports for Fast Reponse and GW
   followup analyses. A lot of code has been taken from Fast Reponse
   Analysis. 

   Author:Sam Hori
   Date:   November 17th, 2026
'''

import numpy as np
import pandas as pd
import os, json, datetime
import logging as log
from dotenv import load_dotenv

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
from . import ReportGenerator

load_dotenv()

class PosteriorReportGenerator(ReportGenerator):
    r"""
    Generate a posterior report from a FastResponseAnalysis object
    """
    _figure_list = []

    def get_report_source(self):
        """ Find base path of the general report TeX document
        Looks for the directory location fast_response, returns 
        fast_response/latex/report_general.tex
        """
        base = os.path.dirname(fast_response.__file__)
        return os.path.join(base, 'latex/posterior_report_general.tex')

    def generate_report(self):
        #just to ensure call to generate report don't create an error
        self.generate_report()
        
    def generate_posterior_report(self):
        """Generate the report! This function does most of the heavy lifting.
        Starts making the source LaTeX file (named r.tex), queries runs,
        makes detector status plots, includes all plots for the analysis,
        makes all tables for analysis parameters, runs, events, etc.
        """

        report_fname = os.path.join(self.dirname, "r_posterior.tex")

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


        s = self.source

        dataset = Datasets[self.analysis._dataset]
        # Make rate plots and put them in analysis directory

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
            print(self.analysisid)
            for plot_name, plot_path in [
                                         ('PriorSkymap', "PriorSkymap.png"),
                                         ('PosteriorSkymap', "PosteriorMap.png"),
                                         ('PosteriorSkymapZoom',"PosteriorMapZoomed.png"),
                                         ('NsProfilePlot',"nsProfilePlot.png"),
                                         ("PosteriorSkymapZoomWNeutrinos","PosteriorMapZoomedWithEvents.png")
                                         ]:
                if os.path.isfile(self.dirname + '/' + plot_path):
                    f.write(
                        r"\newcommand{"+"\\"+plot_name+"}{"+ self.dirname+"/"+
                        plot_path +
                        "}\n"
                    )
    

    def make_pdf(self):
        """Makes the PDF compiled version of the LaTeX report
        using pdflatex, and saves to the analysis folder
        """
        # symlink main report tex file
        reportfname = self.analysisid + "_posterior_report.tex"
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
             self.dirname + '/' + self.analysisid + "_posterior_report.tex"],
            #cwd=self.dirname,
            env = env,
        )

class FastResponsePosteriorReport(PosteriorReportGenerator):
    """Generates a general FastResponse Report for a point source.
    Includes additional plots of the 1D LLH ns scan
    """
    _figure_list = [('nsscan', 'llh_ns_scan.png')]

    def __init__(self, analysis):
        super().__init__(analysis)

    def generate_posterior_report(self):
        super().generate_posterior_report()


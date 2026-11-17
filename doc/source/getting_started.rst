Getting Started
==================


Installation
---------------------
To use the Fast Response Analysis, first clone it from GitHub:

.. code-block:: console

   $ git clone git@github.com:icecube/FastResponseAnalysis.git

We recommend using a virtual environment to install the FRA and its related dependancies. The required libraries are listed in the file requirements.txt, included in the repository. The current iteration of Fast Response uses Python 3.7.5.

.. code-block:: console

   $ python3 -m venv fra_env
   $ source fra_env/bin/activate
   $ pip install -r [PATH_TO_FRA]/requirements.txt

Fast Response can also be pip installed into the venv:

.. code-block:: console

   $ pip install -e .

Fast Response also has dependencies on [private] IceCube software, namely the packages `IceTray <https://github.com/icecube/icetray>`_, `realtime <https://github.com/icecube/realtime>`_ (a parasitic metaproject built onto IceTray) and `Skylab <https://github.com/icecube/skylab>`_, for querying data and analysis tools.

Environment options
^^^^^^^^^^^^^^^^^^^^^
To specify where the output should go, you can set the environment variable: 

.. code-block:: console

   $ export FAST_RESPONSE_OUTPUT=[PATH_TO_OUTPUT_DIR] 

You can also specify where the scripts directory is, using another environment variable:

.. code-block:: console

   $ export FAST_RESPONSE_SCRIPTS=[PATH_TO_SCRIPTS]

Tutorials
---------------------

Externally Triggered Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to perform a short timescale followup using the realtime GFU stream, you need only know:

#. The time window of interest
#. The location of the source (either RA, Dec or localization probability skymap)

In order to run an analysis, navigate to the scripts directory. All you need to do is:

.. code-block:: python

   python run_external_followup.py --name="Fast Response Test" --start="2020-01-01 12:00:00.00" --stop="2020-01-02 12:00:00.0" --ra=250.0 --dec=45.0

Including test (any case) in the name will run with scrambled data, rather than real (UNBLINDED) data. (There is a skymap arg here, that can be passed instead of an RA, Decl - but not all functions are available for skymaps in realtime, so the point source version should be strongly preferred).

After a few seconds you should see the following screen:

.. code-block:: console

   ********************************************************************************
    _____         _     ____                                      
   |  ___|_ _ ___| |_  |  _ \ ___  ___ _ __   ___  _ __  ___  ___ 
   | |_ / _` / __| __| | |_) / _ \/ __| '_ \ / _ \| '_ \/ __|/ _ \
   |  _| (_| \__ \ |_  |  _ <  __/\__ \ |_) | (_) | | | \__ \  __/
   |_|  \__,_|___/\__| |_| \_\___||___/ .__/ \___/|_| |_|___/\___|
                                      |_|                         
    _____     _ _                           
   |  ___|__ | | | _____      ___   _ _ __  
   | |_ / _ \| | |/ _ \ \ /\ / / | | | '_ \ 
   |  _| (_) | | | (_) \ V  V /| |_| | |_) |
   |_|  \___/|_|_|\___/ \_/\_/  \__,_| .__/ 
                                     |_|    
  
   ********************************************************************************
 
   Working on unscrambled (UNBLINDED) data
   ********************************************************************************
   *                                                                              *
   *                            Fast Response Example                             *
   *                                                                              *
   ********************************************************************************
             2020-01-01 12:00:00.000          2020-01-02 12:00:00.000
             RA, dec.: 250.00, 45.000          Extension: 0.0

This printout is meant as a double check for the user, to ensure that the proper source details were passed. This includes a warning that the analysis is being performed on UNBLINDED data. If you want to run an analysis on scrambled data for verification, include the string test anywhere in the name argument of the script.

It will then take about one minute to initialize all relevant objects for the likelihood analysis. If everything ran properly, after messages about time to initialize, you will see a printout with the best-fit ns value as well as the test statistic. If the best-fit ns is zero, then no background trials will be run (p=1.0), otherwise, 1000 background trials will be run to quantify the significance.

After this, an upper limit will be calculated. You should see a printout like:

.. code-block:: console

   Beginning upper limit calculation
   Initializing Point Source Injector
   Found upper limit of 2.44 events

If the analysis was run on a short time window and the best-fit ns was 0, a nice sanity check is that this number should be within about 10% from 2.3 events. Once the analysis is complete, there will be some performance plots that are generated and saved to a folder. A pdf report is also generated in this folder.

Additionally, the following arguments can be passed at the command line:

* extension source extension in degrees for following up extended sources
* skip-events=run_id:event_id for removing events from the analysis (ie when following up high energy events)
* ntrials to change the number of background trials performed in the case of a non-zero TS
* n_per_sig to change the number of signal trials performed when calculating the upper limit

Alert Followup
^^^^^^^^^^^^^^
To run the analysis to follow up an alert event, you will either run run_track_followup.py or run_cascade_followup.py depending on which alert type you are following up. You need to wait until at least 24 hours after the alert event time in order to run the analysis.

In order to run this analysis, you should find the alert event time (sent in the GCN notice, and also written in a GCN Circular, if sent). You'll need to convert the ISO UTC time to MJD (I recommend `this tool <https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xTime/xTime.pl>`_). You'll also need the alert run number and event number (also sent in the GCN Notice).

Once you have these details, you can run the analysis with

.. code-block:: console

   python run_track_followup.py --skymap=/home/followup/output_plots/run{RUNID}.evt{EVENTID}.HESE.skymap_nside_512.fits.gz --time={ALERT_MJD} --gcn_notice_num={GCN_CIRCULAR_NUMBER} --alert_id={RUNID}:{EVENTID}

This will run two analyses, one with a time window of 1000 s and one with a time window of 2 days, both centered on the alert time. It will also remove the alert event from the sample, and it will assume a spectral index of -2.5, which is different than the -2 used for the run_external_followup.py script.

GW Followup
^^^^^^^^^^^^
To run the analysis to follow up a graviational wave event, you will need to navigate to the directory including run_gw_followup.py. You need to know:

#. The name of the GW event
#. The time of the GW
#. A link to the skymap for the GW (can be a url or a path)

As an example, there is a sample map included here: fast_response/sample_skymaps/S191216ap_update.xml. To run this event, the input would look like:

.. code-block:: console

   python run_gw_followup.py --name="S191216ap Update" --time=58833.893 --skymap="https://gracedb.ligo.org/api/superevents/S191216ap/files/LALInference.fits.gz,0"

Additionally, you can choose to pass the argument --allow_neg_ts (bool) if you want to use the convention where a negative TS is allowed. The default is to use the convention TS>=0.

When running the code, you will see a printout like (similar to the externally triggered follow up):

.. code-block:: console

   ********************************************************************************
     ______        __  _____     _ _                           
    / ___\ \      / / |  ___|__ | | | _____      ___   _ _ __  
   | |  _ \ \ /\ / /  | |_ / _ \| | |/ _ \ \ /\ / / | | | '_ \ 
   | |_| | \ V  V /   |  _| (_) | | | (_) \ V  V /| |_| | |_) |
    \____|  \_/\_/    |_|  \___/|_|_|\___/ \_/\_/  \__,_| .__/ 
                                                        |_|    
   
   ********************************************************************************
   Working on unscrambled (UNBLINDED) data
   Grabbing data

Performance plots will be saved and a report with the results generated in the same way as the above analyses. The output directory is specified in the same way as above, with the FAST_RESPONSE_OUTPUT environment variable.

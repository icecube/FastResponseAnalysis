# Fast Response Analysis

## Overview 
The fast response analysis is a way of using established transient likelihood methods to quickly respond to astronomical events of interest in realtime. 

The current implementation houses the code for responding to 
* Generic externally triggered transients (eg. bright GRBs or blazar flares)
* IceCube alert events
* Gravitational Wave events (R. Hussain's analysis)

## Dependencies
To install all dependencies external from IceCube software, we recommend creating a virtual environment

In order to grab events from the i3live database, you will also need a the `realtime` metaproject. `realtime` is a parasitic metaproject, meaning you need a built version of `icetray` before building the `realtime` environment. If you do not have `icetray` installed, see [the DOCS](https://docs.icecube.aq/icetray/main/)

After you have a version of `icetray` built, you will want to do a parasitic build of the realtime project, building off of this version of `icetray`. To do this, navigate to a new directory where you want your realtime project to live, and run

```console
git clone https://github.com/icecube/realtime.git src
mkdir build
cd build
cmake ../src/ -DMETAPROJECT=/path/to/icerec/build/ -DCMAKE_INSTALL_PREFIX=icerec-plus.${OS_ARCH}
make
```

You can now load the realtime project with 
```console
/path/to/realtime/build/env-shell.sh
```

You will also need a version of `skylab` to run this analysis. See [the Skylab GitHub](https://github.com/icecube/skylab) for instructioins on how to install this code and all of its dependencies.

Once you are in your realtime project, you will need to install the relevant dependencies this project requires. We recommend making a virtual environment and installing the relevant dependencies by running the following lines

```console
python3 -m venv fra_env
source fra_env/bin/activate
pip install -r /path/to/fast-response/requirements.txt
```

This will create a virtual environment named `fra_env`, and the `source fra_env/bin/activate` line will activate the environment.

You should then install `fast_response` into this environment by running

```console
pip install fast_response
```

In the future, you will not need to jump through these hoops, and you can load the environment with these lines:

```console
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/setup.sh`
/path/to/realtime/build/env-shell.sh
source /path/to/fra_env/bin/activate
```


## Tutorial — externally triggered analysis
In order to perform a short timescale followup using the realtime GFU stream, you need only know:
1. The time window of interest
2. The location of the source (either RA, Dec or localization probability skymap)

In order to run an analysis, navigate to the directory containing `run_external_followup.py`. All you need to do is
```console
python run_external_followup.py --name="Fast Response Example" --start="2020-01-01 12:00:00.00" --stop="2020-01-02 12:00:00.0" --ra=250.0 --dec=45.0
```
You can either pass the location as an RA, DEC or pass a skymap with `--skymap=/path/to/skymap` (this can be a string with a url or a path if running on the cobalts). After a few seconds, you should see the following screen:

```
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
```
This printout is meant as a double check for the user, to ensure that the proper source details were passed. This includes a warning that the analysis is being performed on UNBLINDED data. If you want to run an analysis on scrambled data for verification, include the string `test` anywhere in the `name` argument of the `run_external_followup.py` script.

It will then take about one minute to initialize all relevant objects for the likelihood analysis. If everything ran properly, after messages about time to initialize, you will see a printout with the best-fit `ns` value as well as the test statistic. If the best-fit `ns` is zero, then no background trials will be run (`p=1.0`), otherwise, 1000 background trials will be run to quantify the significance. 

After this, an upper limit will be calculated. You should see a printout like:
```
Beginning upper limit calculation
Initializing Point Source Injector
Found upper limit of 2.44 events
``` 

If the analysis was run on a short time window and the best-fit `ns` was 0, a nice sanity check is that this number should be within about 10% from 2.3 events. Once the analysis is complete, there will be some performance plots that are generated and saved to a folder. A pdf report is also generated in this folder.

### Options
To specify where the output should go, you can set the environment variable with 
```console
export FAST_RESPONSE_OUTPUT=/path/to/output/
```

Additionally, the following arguments can be passed at the command line:
* `--extension` source extension in degrees for following up extended sources
* `--skip-events=run_id:event_id` for removing events from the analysis (ie when following up high energy events)
* `--ntrials` to change the number of background trials performed in the case of a non-zero TS
* `--document` this is a flag, if present then the official documentation page will be updated upon completion of the analysis

## Tutorial - Alert followup
To run the analysis to follow up an alert event, you will either run `run_track_followup.py` or `run_cascade_followup.py` depending on which alert type you are following up. You need to wait until at least 24 hours after the alert event time in order to run the analysis.

In order to run this analysis, you should find the alert event time (this is shared in the GCN notice which is in the Bacodine emails, and it will also be written in the circular). You'll need to convert the ISO UTC time to MJD (I recommend [this tool](https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xTime/xTime.pl) if you don't have a way of doing it). You'll also need the alert run number and event number. 

Once you have these details, you can run the analysis with

```console
python run_track_followup.py --skymap=/home/followup/output_plots/run{RUNID}.evt{EVENTID}.HESE.skymap_nside_512.fits.gz --time={ALERT_MJD} --gcn_notice_num={GCN_CIRCULAR_NUMBER} --alert_id={RUNID}:{EVENTID}
```

This will run two analyses, one with a time window of 1000 s and one with a time window of 2 days, both centered on the alert time. It will also remove the alert event from the sample, and it will assume a spectral index of -2.5, which is different than the -2 used for the `run_external_followup.py` script.

## Tutorial — Gravitational wave followup
To run the analysis to follow up a graviational wave event, you will need to navigate to the directory including `run_gw_followup.py`. You need to know:
1. The name of the GW event
2. The time of the GW
3. A link to the skymap for the GW (can be a string url link or a path if run on the cobalts

As an example, there is a sample map included here: `fast_response/sample_skymaps/S191216ap_update.xml`. To run this event, the input would look like:
```console
python run_gw_followup.py --name="S191216ap Update" --time=58833.893 --skymap="../sample_skymaps/S191216ap_update.xml"
```
Additionally, you can choose to pass the argument `--allow_neg_ts` (bool) if you want to use the convention where a negative TS is allowed. The default is to use the convention TS>=0.

When running the code, you will see a printout like (similar to the externally triggered follow up):
```
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
```

Performance plots will be saved and a report with the results generated in the same way as the above analyses. The output directory is specified in the same way as above, with the `FAST_RESPONSE_OUTPUT` environment variable.

## Contacts
Please send any questions or feature requests to:
* Jessie Thwaites (thwaites@wisc.edu)
* Alex Pizzuto (apizzuto@icecube.wisc.edu)
* Raamis Hussain (raamis.hussain@icecube.wisc.edu)
* Justin Vandenbroucke (justin.vandenbroucke@wisc.edu)

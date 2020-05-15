# Fast Response Analysis

## Overview 
The fast response analysis is a way of using established transient likelihood methods to quickly respond to astronomical events of interest in realtime. Specifically, the fast response analysis has been used to respond to blazar flares that are discussed on the Astronomer's Telegram, especially bright GRBs (such as GRB190114C), 

## Dependencies
This code relies heavily on the `skylab` analysis software framework, as well as on many modules in the scientific python suite, including `numpy`, `scipy`, `healpy`, `astropy`, and many others. For convenience, a `requirements.txt` file is provided. Right now, only python2 is supported (we hope to add python3 compatibility in the future). You may need to specify this when initializing a virtual environment.

In order to create a virtual environment to run this analysis, I recommend using `virtualenv`, which can be installed via `pip`:

```console
pip install --user virtualenv
```

Then, create a virtual environment with the relevant packages by running
```console
python -m virtualenv fra_env
source fra_env/bin/activate
pip install -r requirements.txt
```

This will create a virtual environment names `fra_env`, and the `source fra_env/bin/activate` line will activate the environment.

## Tutorial
In order to perform a short timescale followup using the realtime GFU stream, you need only know:
1. The time window of interest
2. The location of the source (either RA, Dec or localization probability skymap)

In order to run an analysis, navigate to the directory containing `run.py`. All you need to do is
```console
python run.py --name="Fast Response Example" --start="2020-01-01 12:00:00.00" --stop="2020-01-02 12:00:00.0" --loc="250.0, -45.0"
```
where `loc` is the location of the object, passed either in the format `RA, DEC` in degrees, or to the location of a skymap (this can be a string with a url or a path if running on the cobalts). After a few seconds, you should see the following screen:

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
          RA, dec.: 250.00, -45.000          Extension: 0.0
```
This printout is meant as a double check for the user, to ensure that the proper source details were passed. This includes a warning that the analysis is being performed on UNBLINDED data. If you want to run an analysis on scrambled data for verification, include the string `test` anywhere in the `name` argument of the `run.py` script.

It will then take about one minute to initialize all relevant objects for the likelihood analysis. If everything ran properly, after messages about time to initialize, you will see a printout with the best-fit `ns` value as well as the test statistic. If the best-fit `ns` is zero, then no background trials will be run (`p=1.0`), otherwise, 1000 background trials will be run to quantify the significance. 

After this, an upper limit will be calculated. You should see a printout like:
```
Beginning upper limit calculation
Initializing Point Source Injector
Found upper limit of 2.44 events
``` 

If the analysis was run on a short time window and the best-fit `ns` was 0, a nice sanity check is that this number should be within about 10% from 2.3 events. Once the analysis is complete, there will be some performance plots that are generated and saved to a folder. A pdf report is also generated in this folder, and if you have permissions, it will be copied to a central location and automatically update the [documentation page](https://icecube.wisc.edu/~apizzuto/FastResponse/webpage/). 

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

## Contacts
Please send any questions or feature requests to:
* Alex Pizzuto (apizzuto@icecube.wisc.edu)
* Raamis Hussain (raamis.hussain@icecube.wisc.edu)
* Justin Vandenbroucke (justin.vandenbroucke@wisc.edu)

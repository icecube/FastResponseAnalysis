# Fast Response Analysis

## Overview 
The fast response analysis is a way of using established transient likelihood methods to quickly respond to astronomical events of interest in realtime. Specifically, the fast response analysis has been used to respond to blazar flares that are discussed on the Astronomer's Telegram, especially bright GRBs (such as GRB190114C), 

## Dependencies
This code relies heavily on the `skylab` analysis software framework, as well as on many modules in the scientific python suite, including `numpy`, `scipy`, `healpy`, `astropy`, and many others. For convenience, a `requirements.txt` file is provided for each stable release. This release was written and tested with python3 (`py3-v4.1.0` on the cobalts). 

In order to grab events from the i3live database, you will also need a the `realtime` metaproject. The lines below walkthrough how to create a virtual environment with all of these dependencies, assuming the user is working on the IceCube filesystem. First, navigate to the location where you want to put icerec, and run

### Note on Python Versions
The trunk (and any release including and after v02-00) use python3. v01-00 uses python2, and as such, if you would like to use that version, you may have to adjust the installation below to point to relevant python2 environments.

```console
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/setup.sh`
svn co http://code.icecube.wisc.edu/svn/meta-projects/combo/releases/V01-00-00/ src
mkdir build 
cd build
cmake ../src
make
```

After you have a version of icerec built, you will want to do a parasitic build of the realtime project, building off of this version of icerec. To do this, navigate to a new directory where you want your realtime project to live, and run

```console
svn co http://code.icecube.wisc.edu/svn/meta-projects/realtime/trunk/ src
mkdir build 
cd build
cmake ../src/ -DMETAPROJECT=/path/to/icerec/build/ -DCMAKE_INSTALL_PREFIX=combo-plus.${OS_ARCH}
make
```

You can now load the realtime project with 
```console
/path/to/realtime/build/env-shell.sh
```

Once you are in your realtime project, you will need to install the relevant dependencies this project requires. We recommend making a virtual environment and installing the relevant dependencies by running the following lines

```console
python3 -m venv fra_env
source fra_env/bin/activate
pip install -r /path/to/fast-response/requirements.txt
```

This will create a virtual environment names `fra_env`, and the `source fra_env/bin/activate` line will activate the environment.

In the future, you will not need to jump through these hoops, and you can load the environment with these lines:

```console
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/setup.sh`
/path/to/realtime/build/env-shell.sh
source /path/to/fra_env/bin/activate
```

If you would like to be able to import fast response tools from any directory in the future, you must append to your python path with 
```console
export PYTHONPATH=$PYTHONPATH:/path/to/fast-response
```

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

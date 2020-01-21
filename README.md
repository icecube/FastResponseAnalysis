# Fast Response Analysis

## Overview 
The fast response analysis is a way of using established transient likelihood methods to quickly respond to astronomical events of interest in realtime. Specifically, the fast response analysis has been used to respond to blazar flares that are discussed on the Astronomer's Telegram, especially bright GRBs (such as GRB190114C), 

## Dependencies
This code relies heavily on the `skylab` analysis software framework, as well as on many modules in the scientific python suite, including `numpy`, `scipy`, `healpy`, `astropy`, and many others. For convenience, a `requirements.txt` file is provided.

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

It will then take about one minute to initialize all relevant objects for the likelihood analysis. If everything ran properly, after messages about time to initialize, you will see 

```
TS = 0.0
ns = 0.0

```


## Code


## Contacts
* Alex Pizzuto (apizzuto@icecube.wisc.edu)
* Raamis Hussain ()
* Justin Vandenbroucke ()
* Kevin Meagher ()

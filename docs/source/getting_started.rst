Getting Started
==================

[Under construction]

Installing FRA
---------------------


Tutorials
---------------------
In order to perform a short timescale followup using the realtime GFU stream, you need only know:

# The time window of interest
# The location of the source (either RA, Dec or localization probability skymap)

In order to run an analysis, navigate to the scripts directory. All you need to do is:
.. code-block:: python

   python run_external_followup.py --name="Fast Response Test" --start="2020-01-01 12:00:00.00" --stop="2020-01-02 12:00:00.0" --ra=250.0 --dec=45.0

Including test (any case) in the name will run with scrambled data, rather than unblinded data.

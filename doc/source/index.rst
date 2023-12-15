.. FastResponseAnalysis documentation master file, created by
   sphinx-quickstart on Mon Dec  4 10:43:46 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

IceCube Fast Response Analysis (FRA)
=====================================

::
   
   ********************************************************************************
    _____         _     ____                                      
   |  ___|_ _ ___| |_  |  _ \ ___  ___ _ __   ___  _ __  ___  ___ 
   | |_ / _` / __| __| | |_) / _ \/ __| '_ \ / _ \| '_ \/ __|/ _ \
   |  _| (_| \__ \ |_  |  _ <  __/\__ \ |_) | (_) | | | \__ \  __/
   |_|  \__,_|___/\__| |_| \_\___||___/ .__/ \___/|_| |_|___/\___|
                                      |_|                         
       _                _           _     
      / \   _ __   __ _| |_   _ ___(_)___ 
     / _ \ | '_ \ / _` | | | | / __| / __|
    / ___ \| | | | (_| | | |_| \__ \ \__ \
   /_/   \_\_| |_|\__,_|_|\__, |___/_|___/
                          |___/           

   ********************************************************************************

The **Fast Response Analysis** is a tool for rapid response to astrophysical transients with the IceCube Neutrino Observatory. It uses a transient maximum likelihood method to rapidly respond to astrophysical events in real time.

The current implementation houses the code for responding to:

* Generic externally triggered transients (eg. bright GRBs, AGN flares, novae, etc)
* IceCube alert events
* Gravitational Wave events

This documentation describes the Fast Response Package. For a description of the physics, see `R. Abbasi et al 2021 ApJ 910 4 <https://iopscience.iop.org/article/10.3847/1538-4357/abe123>`_ (for external transients), `R. Abbasi et al 2023 ApJ 951 45 <https://iopscience.iop.org/article/10.3847/1538-4357/acd2ca>`_ (for alert events), and `R. Abbasi et al 2023 ApJ 944 80 <https://iopscience.iop.org/article/10.3847/1538-4357/aca5fc>`_ (for gravitational wave events).

(Top graphic made with `PyFiglet <https://pypi.org/project/pyfiglet/>`_)

Contents
============
.. toctree::
   :maxdepth: 2

   getting_started
   fra
   alert
   ext
   gw
   reports
   web_utils
   plotting

Index
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

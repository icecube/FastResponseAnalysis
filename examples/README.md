# Making paper plots

The scripts in this directory will generate all of the plots that are included in the fast response analysis paper (except for zoomed skymaps, for which you should run the fast response analysis with analysis details for each as given by the entries in https://roc.icecube.wisc.edu/internal/fast_response/).

There are three scripts that make paper plots. They are:
* `analysis_paper_plots.py`: Generates all sensitivity and median significance plots
* `followup_paper_plots.py`: Generates SEDs for GRB190114C and PKS0346-27
* `FRB_population_limits.py`: Integrates over a standard candle luminosity distribution of FRBs to calculate an overall limit on SGR1935+2154-like FRBs (requires a working version of [flarestack](https://github.com/IceCubeOpenSource/flarestack))
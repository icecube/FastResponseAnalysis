Plotting
=======================================================
.. note::
   Functions for making plots, both analysis and detector plots, as well as sensitivity calculations. These are used internally in the FRA package. 

Plotting functions
-------------------------------------------------------
.. autofunction:: fast_response.plotting_utils.plot_zoom

.. autofunction:: fast_response.plotting_utils.plot_color_bar

.. autofunction:: fast_response.plotting_utils.plot_labels

.. autofunction:: fast_response.plotting_utils.plot_events

.. autofunction:: fast_response.plotting_utils.load_plotting_settings

.. autofunction:: fast_response.plotting_utils.contour

.. autofunction:: fast_response.plotting_utils.plot_contours


Detector status plots
-------------------------------------------------------
.. autofunction:: fast_response.make_ontime_plots.time_axis
   
.. autofunction:: fast_response.make_ontime_plots.time_series

.. autofunction:: fast_response.make_ontime_plots.make_rate_plots


Sensitivity functions
------------------------
.. autofunction:: fast_response.sensitivity_utils.find_nearest_idx

.. autofunction:: fast_response.sensitivity_utils.find_nearest

.. autofunction:: fast_response.sensitivity_utils.deltaPsi

.. autofunction:: fast_response.sensitivity_utils.deltaPsi2

.. autofunction:: fast_response.sensitivity_utils.sensitivity_fit

.. autofunction:: fast_response.sensitivity_utils.erfunc

.. autofunction:: fast_response.sensitivity_utils.chi2cdf

.. autofunction:: fast_response.sensitivity_utils.incomplete_gamma

.. autofunction:: fast_response.sensitivity_utils.poissoncdf

.. autofunction:: fast_response.sensitivity_utils.binomial_error

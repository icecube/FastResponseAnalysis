AlertFollowup
=======================================================
.. note::
   Class for followup to Alert events, using a spatial prior and a fixed spectral index of 2.5 in the likelihood.

General Alert followup class
------------------------------
.. autoclass:: fast_response.AlertFollowup.AlertFollowup
   :members: generate_report, ps_sens_range, run_background_trials, sens_range_plot, write_circular

Track Followups (Gold and Bronze)
----------------------------------
.. autoclass:: fast_response.AlertFollowup.TrackFollowup
   :members:

Cascade Followups
-------------------------
.. autoclass:: fast_response.AlertFollowup.CascadeFollowup
   :members:

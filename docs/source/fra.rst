FastResponseAnalysis
=======================================================
.. note::
   Shared base class, and classes for both skymap and point-source type followups, with shared functions for each.

Base FastResponseAnalysis class
--------------------------------
.. autoclass:: fast_response.FastResponseAnalysis.FastResponseAnalysis
   :members:

Base PriorFollowup class
---------------------------
.. autoclass:: fast_response.FastResponseAnalysis.PriorFollowup
   :members: dec_skymap_range, find_coincident_events, format_skymap, initialize_injector, ipixs_in_percentage, make_dNdE, run_background_trials, unblind_TS

Base PointSourceFollowup class
---------------------------------
.. autoclass:: fast_response.FastResponseAnalysis.PointSourceFollowup
   :members:

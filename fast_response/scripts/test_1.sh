python test_posterior_cascade.py --skymap=https://roc.icecube.wisc.edu/public/hese_cascades/hese_60505_run00139644.evt000063900267.fits --time=60505.623363 --alert_id=139644:63900267
python test_posterior_cascade.py --skymap=https://roc.icecube.wisc.edu/public/hese_cascades/hese_60505_run00139644.evt000063900267.fits --time=60505.623363 --alert_id=139644:63900267 --prior=jeffries


python test_posterior_code.py --skymap=https://gracedb.ligo.org/apiweb/superevents/S240716b/files/Bilby.fits.gz --time=60000.00 \
    --name="Posterior_Test_2inj" --tw=1000 --n_inj=2

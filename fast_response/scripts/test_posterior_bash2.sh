#python test_posterior_code.py --skymap=https://gracedb.ligo.org/apiweb/superevents/S240716b/files/Bilby.fits.gz --time=60000.00 --name="PosteriorTest_RealEventRandomDate_1000s_2" --tw=1000
#https://gracedb.ligo.org/superevents/S240716b/files/
#astropy.time.Time(1384415586.09,format='gps').mjd

python test_posterior_code.py --skymap=https://gracedb.ligo.org/apiweb/superevents/S240716b/files/Bilby.fits.gz --time=60000.00 \
    --name="Posterior_Test_0inj" --tw=1000 --n_inj=0

python test_posterior_code.py --skymap=https://gracedb.ligo.org/apiweb/superevents/S240716b/files/Bilby.fits.gz --time=60000.00 \
    --name="Posterior_Test_1inj" --tw=1000 --n_inj=1
    
python test_posterior_code.py --skymap=https://gracedb.ligo.org/apiweb/superevents/S240716b/files/Bilby.fits.gz --time=60000.00 \
    --name="Posterior_Test_2inj" --tw=1000 --n_inj=2

python test_posterior_code.py --skymap=https://gracedb.ligo.org/apiweb/superevents/S240716b/files/Bilby.fits.gz --time=60000.00 \
    --name="Posterior_Test_3inj" --tw=1000 --n_inj=3

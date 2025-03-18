#python test_posterior_code.py --skymap=https://gracedb.ligo.org/apiweb/superevents/S240716b/files/Bilby.fits.gz --time=60000.00 --name="PosteriorTest_RealEventRandomDate_1000s_2" --tw=1000
#https://gracedb.ligo.org/superevents/S240716b/files/
#astropy.time.Time(1384415586.09,format='gps').mjd

#S231119u-4-Update
python test_posterior_code.py --skymap=https://gracedb.ligo.org/api/superevents/S231119u/files/bayestar.fits.gz \
                            --time=60267.328334432874 --name="S231119u_Test_Posterior" --tw=1000
#S230928cb-4-Update
python test_posterior_code.py --skymap=https://gracedb.ligo.org/api/superevents/S230928cb/files/bayestar.fits.gz \
                            --time=60215.915591805555 --name="S230928cb_Test_Posterior" --tw=1000
#S230518h-4-Update
python test_posterior_code.py --skymap=https://gracedb.ligo.org/api/superevents/S230518h/files/bayestar.fits.gz \
                            --time=60082.54106674768 --name="S230518h_Test_Posterior" --tw=1000
#S230522a-2-Preliminary
python test_posterior_code.py --skymap=https://gracedb.ligo.org/api/superevents/S230522a/files/bayestar.fits.gz \
                            --time=60086.40144832176 --name="S230522a_Test_Posterior" --tw=1000

#S230529ay-4-Update
python test_posterior_code.py --skymap=https://gracedb.ligo.org/api/superevents/S230529ay/files/bayestar.fits.gz \
                            --time=60093.76042530093 --name="S230529ay_Test_Posterior" --tw=1000
#S240615dg-4-Update
python test_posterior_code.py --skymap=https://gracedb.ligo.org/api/superevents/S240615dg/files/bayestar.fits.gz \
                            --time=60476.48357319445 --name="S240615dg_Test_Posterior" --tw=1000
#S231029y-4-Update
python test_posterior_code.py --skymap=https://gracedb.ligo.org/api/superevents/S231029y/files/bayestar.fits.gz \
                            --time=60246.46885131944 --name="S231029y_Test_Posterior" --tw=1000
#S230904n-4-Update
python test_posterior_code.py --skymap=https://gracedb.ligo.org/api/superevents/S230904n/files/bayestar.fits.gz \
                            --time=60191.21542972222 --name="S230904n_Test_Posterior" --tw=1000

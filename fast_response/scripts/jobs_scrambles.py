#!/usr/bin/env python
import itertools
import subprocess
import processing_tools as tools
from glob import glob

new_examples = [
    'IGWN-GWTC5p0-29ebe06b7_25-GW240531_075248-NRSur7dq4_Skymap_PEDataRelease.converted.fits',
    'IGWN-GWTC5p0-29ebe06b7_25-GW240910_103535-IMRPhenomXPNR_Skymap_PEDataRelease.converted.fits',
    'IGWN-GWTC5p0-29ebe06b7_25-GW250108_152221-IMRPhenomXPHM_SpinTaylor_Skymap_PEDataRelease.converted.fits',
    'IGWN-GWTC5p0-29ebe06b7_25-GW240908_125134-IMRPhenomXPHM_SpinTaylor_Skymap_PEDataRelease.converted.fits',
    'IGWN-GWTC5p0-29ebe06b7_25-GW241009_220455-IMRPhenomXPNR_Skymap_PEDataRelease.converted.fits',
    'IGWN-GWTC5p0-29ebe06b7_25-GW241002_030559-IMRPhenomXPNR_Skymap_PEDataRelease.converted.fits',
    'IGWN-GWTC5p0-29ebe06b7_25-GW240615_113620-SEOBNRv5PHM_Skymap_PEDataRelease.converted.fits',
    'IGWN-GWTC5p0-29ebe06b7_25-GW250114_082203-SEOBNRv5PHM_Skymap_PEDataRelease.converted.fits',
    'IGWN-GWTC5p0-29ebe06b7_25-GW240703_191355-IMRPhenomXPHM_SpinTaylor_Skymap_PEDataRelease.converted.fits',
    'IGWN-GWTC5p0-29ebe06b7_25-GW240930_035959-SEOBNRv5PHM_Skymap_PEDataRelease.converted.fits',
    'IGWN-GWTC5p0-29ebe06b7_25-GW240908_082628-IMRPhenomXPHM_SpinTaylor_Skymap_PEDataRelease.converted.fits',
    'IGWN-GWTC5p0-29ebe06b7_25-GW241129_021832-SEOBNRv5PHM_Skymap_PEDataRelease.converted.fits',
    ]
gwtc_examples = []
for gw_name in new_examples:
    gwtc_examples.append(glob(f'/data/user/chraab/GWTC-5.0/skymaps/*{gw_name}*')[0])
all_gwtc = glob('/data/user/chraab/GWTC-5.0/skymaps/*converted.fits')
config = dict(
    script='/data/user/chraab/metaprojects/realtime/FastResponseAnalysis/fast_response/scripts/dnnc_scrambles.py',
    duration=1000./86400,
    #bkg=[0.21],
    #nside=128,
    skymap=gwtc_examples,
    #seed=0,
    seed=list(range(1000, 6000, 1000)),
    nseed=100,
    index=2,
    outdir='/data/user/chraab/fast_response/multisample/extended_archival',
    start = '2026-06-30',
)
if 'fix_index' not in config:
    config['outdir'] += '_floating'


with open('jobs.txt', 'w') as f:
    for _config in tools.generate_arg_combinations(config):
        script = _config.pop('script')
        arg_string = tools.build_arg_string(_config)
        f.write(' '.join([script, arg_string, '\n']))

subprocess.run('condor_submit job_alma9.sub', shell=True, check=True)

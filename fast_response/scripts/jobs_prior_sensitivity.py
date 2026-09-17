#!/usr/bin/env python
from glob import glob
import pandas as pd
import itertools
import subprocess
import numpy as np
import processing_tools as tools

examples = ['GW240514_121713-NRSur7dq4_Skymap_PEDataRelease.converted.fits',
 'GW240924_000316-IMRPhenomXPHM_SpinTaylor_Skymap_PEDataRelease.converted.fits',
 'GW240527_230910-SEOBNRv5PHM_Skymap_PEDataRelease.converted.fits',
 'GW240915_001357-IMRPhenomXPHM_SpinTaylor_Skymap_PEDataRelease.converted.fits',
 'GW240505_133552-NRSur7dq4_Skymap_PEDataRelease.converted.fits',
 'GW240511_031507-SEOBNRv5PHM_Skymap_PEDataRelease.converted.fits',
 'GW240621_214041-NRSur7dq4_Skymap_PEDataRelease.converted.fits',
 'GW240527_183429-IMRPhenomXPNR_Skymap_PEDataRelease.converted.fits',
 'GW240531_040326-SEOBNRv5PHM_Skymap_PEDataRelease.converted.fits',
 'GW240825_055146-IMRPhenomXPNR_Skymap_PEDataRelease.converted.fits']
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
gwtc_all = glob('/data/user/chraab/GWTC-5.0/skymaps/*converted.fits')
pslim = glob('/data/user/sahori/GWAnalysisWDNNCascade_Storage/PointSourceLimitTesting/FakeGWs/PSLim*fits')

config = dict(
    script='/data/user/chraab/metaprojects/realtime/FastResponseAnalysis/fast_response/scripts/dnnc_prior_sensitivity_scan.py',
    index = 2.0,
    duration=[1000./3600/24],
    # fix_index='',
    start = '2025-06-30',
    n_iter=1000,
    # alpha=1.35e-3, beta=0.5, mu=list(np.arange(0.55, 0.95, 0.02)),
    # mu = list(np.arange(0.5, 2.0, 0.1)),
    #alpha=0.5, beta=0.9, mu=list(np.arange(2.2, 3.6, 0.04)),
    alpha=0.5, beta=0.9, mu=list(np.arange(4, 10, 0.5)),
    # mu = list(np.arange(2.1, 3.5, 0.1)),
    #skymap = pslim[:5],
    # skymap = gwtc_examples + pslim[:5],
    skymap = gwtc_examples + pslim,
    out='/data/user/chraab/fast_response/multisample/extended_archival',
)
if "fix_index" not in config:
    config["out"] += "_floating"

with open('jobs.txt', 'w') as f:
    for _config in tools.generate_arg_combinations(config):
        script = _config.pop('script')
        arg_string = tools.build_arg_string(_config)
        f.write(' '.join([script, arg_string, '\n']))

subprocess.run('condor_submit job_alma9.sub', shell=True, check=True)

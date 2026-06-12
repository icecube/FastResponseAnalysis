#!/bin/bash
eval $(/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/setup.sh)
source /data/user/chraab/metaprojects/realtime/venv_py3-v4.4.2/bin/activate
/data/user/chraab/metaprojects/realtime/build_alma9/env-shell.sh python "$@"


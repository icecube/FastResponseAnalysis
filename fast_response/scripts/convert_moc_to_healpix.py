import numpy as np
import matplotlib as mpl
mpl.use('agg')
import mhealpy as mhp
import healpy as hp
from astropy.table import QTable
from astropy import units as u
import argparse

parser = argparse.ArgumentParser(description='GW Followup')
parser.add_argument('--skymap', type=str, default=None,
                    help='path to skymap (should be the *.multiorder.fits downloaded from GraceDB)')
args = parser.parse_args()

# Read skymap
skymap = QTable.read(args.skymap)
# convert to a probability in each pixel by multiplying by pix area
s = [skymap[i]['PROBDENSITY'].to_value(u.deg**-2)*hp.nside2pixarea(mhp.uniq2nside(skymap[i]['UNIQ']),degrees=True)
      for i in range(len(skymap))]

# put into mhealpy & get flattened map
m = mhp.HealpixMap(data=s, uniq=skymap['UNIQ'])
#returns an nside of 256 for new map
new_map = m.rasterize(256, 'NESTED')

# save the new map
new_map.write_map(args.skymap.replace('multiorder','converted'), overwrite=True)
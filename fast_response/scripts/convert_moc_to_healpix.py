#!/usr/bin/env python

import sys
import matplotlib as mpl
mpl.use('agg')
import mhealpy as mhp
import healpy as hp
from astropy.table import QTable
from astropy import units as u
import argparse

parser = argparse.ArgumentParser(description='GW Followup')
parser.add_argument('--skymap', dest="skymap_file", type=str, required=True,
                    help='path to skymap (should be the *.multiorder.fits downloaded from GraceDB)')
parser.add_argument('--nside', type=int, default=512,
                    help='nside to be used with the skymap (default=512)')


def main(skymap_file: str, nside: int=512) -> str:
    """Convert a multiorder map to a flat one

    Parameters
    ----------
    skymap_file : str
        Path of the multiorder map
    nside : int, optional
        Output nside, by default 512

    Returns
    -------
    str
        Output file written, will have "multiorder" replaced by "converted"
        or be identical to `skymap` in absence of that substring.
    """
    # Read skymap
    skymap = QTable.read(skymap_file)
    # convert to a probability in each pixel by multiplying by pix area
    s = [skymap[i]['PROBDENSITY'].to_value(u.deg**-2)*hp.nside2pixarea(mhp.uniq2nside(skymap[i]['UNIQ']),degrees=True) for i in range(len(skymap))]

    # put into mhealpy & get flattened map
    m = mhp.HealpixMap(data=s, uniq=skymap['UNIQ'])
    #returns a new map
    new_map = m.rasterize(nside, 'NESTED')

    # save the new map
    outfile = skymap_file.replace('multiorder','converted')
    new_map.write_map(outfile, overwrite=True)
    return outfile

if __name__ == "__main__":
    args = parser.parse_args()
    sys.exit(main(**vars(args)))
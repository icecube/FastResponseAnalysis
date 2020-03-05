#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import healpy as hp
import pickle
import matplotlib.pyplot as plt
from astropy.io import fits
from glob import glob
from time import time
import meander
import argparse

parser = argparse.ArgumentParser(description='Make Contours from probability maps')
parser.add_argument('--ind', type=int, default=0,
                    help='Index of skymap')
parser.add_argument('--nsp', type=int, default=6,
                    help='nside = 2**nsp')
args = parser.parse_args()

base_path = '/data/user/steinrob/millipede_scan_archive/fits_v3_prob_map/'
files = glob(base_path + '*.fits')
nside = 2**args.nsp
proportions = [0.5, 0.68, 0.90, 0.99]
f = files[args.ind]
name = f[60:-5]

def plot_contours(proportions,samples):
    r''' Plot containment contour around desired level.
    E.g 90% containment of a PDF on a healpix map

    Parameters:
    -----------
    proportions: list
        list of containment level to make contours for.
        E.g [0.68,0.9]
    samples: array
        array of values read in from healpix map
        E.g samples = hp.read_map(file)
    Returns:
    --------
    theta_list: list
        List of arrays containing theta values for desired contours
    phi_list: list
        List of arrays containing phi values for desired contours
    '''

    levels = []
    t0 = time()
    sorted_samples = list(reversed(list(sorted(samples))))
    t1 = time()
    nside = hp.pixelfunc.get_nside(samples)
    sample_points = np.array(hp.pix2ang(nside,np.arange(len(samples)))).T
    for proportion in proportions:
        level_index = (np.cumsum(sorted_samples) > proportion).tolist().index(True)
        level = (sorted_samples[level_index] +
                (sorted_samples[level_index+1] if level_index+1<len(samples) else 0))/2.0
        levels.append(level)
    msk = samples != 0
    samples = samples[msk]
    sample_points = sample_points[msk]
    contours_by_level = meander.spherical_contours(sample_points, samples, levels)

    theta_list = []; phi_list=[]
    for contours in contours_by_level:
        for contour in contours:
            theta, phi = contour.T
            phi[phi<0] += 2.0*np.pi
            theta_list.append(theta)
            phi_list.append(phi)

    return theta_list, phi_list
            
def skymap_from_fits(path, nside = 2**8):
    a = fits.open(path)
    probs = hp.pixelfunc.ud_grade(a[0].data, nside)
    probs = probs/np.sum(probs)
    return probs

probs = skymap_from_fits(f, nside = nside)
#probs = hp.pixelfunc.ud_grade(skymap, nside, power=-2)
print('Making contours')
theta, phi = plot_contours(proportions,probs)
print('Made contours')
contours = {str(level): (th, ph) for level, th, ph in zip(proportions, theta, phi)}

with open('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/contours/contour_{}_nsp_{}.pickle'.format(name, args.nsp), 'wb') as f:
    pickle.dump(contours, f)

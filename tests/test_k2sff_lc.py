""" Test the K2SFF light curve files."""

from __future__ import print_function, division
import os
import pytest

import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits

from .. import source_lc

data = os.path.join(os.path.dirname(__file__),"k2sff_211889983.fits")
# Extract the original file for testing

def test_k2sff_lc():
    lc = source_lc.source.from_k2sff(data,1)

def test_header():
    lc = source_lc.source.from_k2sff(data,0)

def test_other_ext():
    lc = source_lc.source.from_k2sff(data,8)

#"k2sff_211887567.fits",

def test_k2sff_epochs():
    filenames = ["k2sff_211889983.fits",
                 "k2sff_211892153.fits", "k2sff_211892173.fits"]
    nstars = len(filenames)
    star_list = []
    ep_lengths = np.zeros(nstars,int)
    for i in range(nstars):
        star = source_lc.source.from_k2sff(filenames[i])
        star_list.append(star)
        ep_lengths[i] = len(star.epochs)

    assert(np.all(ep_lengths==ep_lengths[0]))

def test_k2sff_epochs_nofix():
    filenames = ["k2sff_211889983.fits",
                 "k2sff_211892153.fits", "k2sff_211892173.fits"]
    nstars = len(filenames)
    star_list = []
    ep_lengths = np.zeros(nstars,int)
    for i in range(nstars):
        star = source_lc.source.from_k2sff(filenames[i])
        star_list.append(star)
        ep_lengths[i] = len(star.epochs)

    new_list = source_lc.fix_epochs(star_list)

def test_fix_epoch_lengths():
    filenames = ["k2sff_211887567.fits","k2sff_211889983.fits",
                 "k2sff_211892153.fits", "k2sff_211892173.fits"]
    nstars = len(filenames)
    star_list = []
    ep_lengths = np.zeros(nstars,int)
    for i in range(nstars):
        star = source_lc.source.from_k2sff(filenames[i])
        star_list.append(star)
        ep_lengths[i] = len(star.epochs)

    if np.all(ep_lengths==ep_lengths[0])==False:
        new_list = source_lc.fix_epochs(star_list)

    for i,star in enumerate(new_list):
        ep_lengths[i] = len(star.epochs)

    assert(np.all(ep_lengths==ep_lengths[0]))

def test_fix_epoch_lengths2():
    filenames = ["k2sff_211887567.fits","k2sff_211889983.fits",
                 "k2sff_211892153.fits", "k2sff_211892173.fits"]
    nstars = len(filenames)
    star_list = []
    ep_lengths = np.zeros(nstars,int)
    for i in range(nstars):
        star = source_lc.source.from_k2sff(filenames[i])
        star_list.append(star)
        ep_lengths[i] = len(star.epochs)

    if np.all(ep_lengths==ep_lengths[0])==False:
        new_list = source_lc.fix_epochs(star_list)

    ep_lengths2 = np.zeros(nstars,int)
    for i,star in enumerate(new_list):
        ep_lengths2[i] = len(star.epochs)

    # None of these should have gotten longer than the original
    # longest array
    assert(np.all(ep_lengths2==max(ep_lengths)))

def test_fix_epochs():
    filenames = ["k2sff_211887567.fits","k2sff_211889983.fits",
                 "k2sff_211892153.fits", "k2sff_211892173.fits"]
    nstars = len(filenames)
    star_list = []
    ep_lengths = np.zeros(nstars,int)
    for i in range(nstars):
        star = source_lc.source.from_k2sff(filenames[i])
        star_list.append(star)
        ep_lengths[i] = len(star.epochs)
    old_stars = star_list

    if np.all(ep_lengths==ep_lengths[0])==False:
        new_list = source_lc.fix_epochs(star_list)

    ep_check = np.zeros(nstars,bool)
    for i,star in enumerate(new_list):
        ostar = old_stars[i]
        check_diff = star.epochs[star.flags==0] - ostar.epochs[ostar.flags==0]
        ep_check[i] = np.all(check_diff<1e-10)

    assert(np.all(ep_check))

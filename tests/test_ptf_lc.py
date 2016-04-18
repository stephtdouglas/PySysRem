""" Test the light curve files."""

from __future__ import print_function, division
import os
import pytest

import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits

from .. import source_lc

data = os.path.join(os.path.dirname(__file__),"lc_00001_calib.txt")
# Extract the original file for testing
lc_orig = np.loadtxt(data,usecols=np.arange(8))

def test_ptf_lc():
    lc = source_lc.source.from_ptf(data)

def test_ptf_fname_str():
    lc = source_lc.source.from_ptf(data)
    assert type(lc.filename) is str

def test_ptf_filename():
    lc = source_lc.source.from_ptf(data)
    assert lc.filename==data

def test_ptf_epochs():
    lc = source_lc.source.from_ptf(data)
    assert np.all(lc.epochs==lc_orig[:,0])

def test_ptf_mags():
    lc = source_lc.source.from_ptf(data)
    assert np.all(lc.mags==lc_orig[:,1])

# def test_ptf_mag_errs():
#     lc = source_lc.source.from_ptf(data)
#     assert np.all(lc.magerrors==lc_orig[:,2])
#
# def test_ptf_sky():
#     lc = source_lc.source.from_ptf(data)
#     assert np.all(lc.xs==lc_orig[:,3])
#     assert np.all(lc.ys==lc_orig[:,4])
#
# def test_ptf_chips():
#     lc = source_lc.source.from_ptf(data)
#     assert np.all(lc.xschip==lc_orig[:,5])
#     assert np.all(lc.yschip==lc_orig[:,6])

def test_ptf_flags():
    lc = source_lc.source.from_ptf(data)
    assert np.all(lc.flags==lc_orig[:,7])

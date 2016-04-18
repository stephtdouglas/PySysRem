""" Test the sysrem file."""

from __future__ import print_function, division
import os
import pytest

import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits

from .. import source_lc
from .. import sysrem

def test_pass():
    pass

fake_epochs = np.linspace(0,10,100)
fake_flux = np.sin(fake_epochs)
fake_flags = np.zeros(100,int)
fake_bad = np.array([3,27,52,88])
fake_flags[fake_bad] = 1
true_period = 2*np.pi
check_res = np.median(fake_flux[~fake_bad])
fake_errors = np.ones_like(fake_flux)
fake_else = np.zeros_like(fake_flags)

base_lc = source_lc.source(fake_flux, fake_errors, fake_epochs, fake_flags,"")

def test_matrix_err():
    star_list = []
    for i in range(3):
        star = source_lc.source(fake_flux, fake_errors, fake_epochs, fake_flags,
                                str(i))
        star_list.append(star)

    res, err, meds, slist = sysrem.generate_matrix(star_list)

    check_err = (err==1) | (err>10**19)

    assert np.all(check_err)

def test_res_shape():
    star_list = []
    for i in range(3):
        star = source_lc.source(fake_flux, fake_errors, fake_epochs, fake_flags,
                                str(i))
        star_list.append(star)

    res, err, meds, slist = sysrem.generate_matrix(star_list)

    check_shape = np.shape(res)

    assert check_shape[0]==3
    assert check_shape[1]==100

def test_array_len():
    star_list = []
    for i in range(3):
        star = source_lc.source(fake_flux+i, fake_errors, fake_epochs, fake_flags,
                                "faketest/star_{0}.txt".format(i))
        star_list.append(star)

    res, err, meds, slist = sysrem.generate_matrix(star_list)

    len_orig = len(slist[0].orig_epochs)
    assert len_orig==100

def test_sysrem_fakes():
    star_list = []
    for i in range(3):
        star = source_lc.source(fake_flux+i, fake_errors, fake_epochs, fake_flags,
                                "faketest/star_{0}.txt".format(i))
        star_list.append(star)

    res, err, meds, slist = sysrem.generate_matrix(star_list)
    sysrem.sysrem(res, err, meds, slist)

    fake_lc = np.loadtxt("faketest/star_1.sysrem.txt",usecols=np.arange(4))
    lenfake = len(fake_lc[:,0])
    # The output light curve should be as long as the entire original lc
    assert lenfake==100

def test_real_stars():
    filenames = ["lc_00001_calib.txt", "lc_00002_calib.txt",
                 "lc_00003_calib.txt", ]
    star_list = []
    for i in range(3):
        star = source_lc.source.from_ptf(filenames[i])
        star_list.append(star)
    residuals, errors, meds, star_list2 = sysrem.generate_matrix(star_list)
    sysrem.sysrem(residuals, errors, meds, star_list2)

    old_lc = np.loadtxt("lc_00001_sysrem.txt",usecols=np.arange(8))
    old_flux = old_lc[:,1]
    new_lc = np.loadtxt("lc_00001_calib.sysrem.txt",usecols=np.arange(4))
    new_flux = new_lc[:,1]

    flux_diff = np.abs(old_flux - new_flux)
    check_diff = (flux_diff<1e-10) | (new_flux<=-9999)

    assert np.all(check_diff)

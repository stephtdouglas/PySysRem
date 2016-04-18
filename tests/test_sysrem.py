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

    res, err, fnames = sysrem.generate_matrix(star_list)

    check_err = (err==1) | (err>10**19)

    assert np.all(check_err)

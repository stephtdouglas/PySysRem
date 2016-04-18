""" Test the K2SFF light curve files."""

from __future__ import print_function, division
import os
import pytest

import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits

from .. import source_lc

data = os.path.join(os.path.dirname(__file__),"k2sff_211748286.fits")
# Extract the original file for testing

def test_k2sff_lc():
    lc = source_lc.source.from_k2sff(data,1)

def test_header():
    lc = source_lc.source.from_k2sff(data,0)

def test_other_ext():
    lc = source_lc.source.from_k2sff(data,8)

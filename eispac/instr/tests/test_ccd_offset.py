import pytest

import astropy.units as u
import numpy as np

from eispac.instr import ccd_offset


@pytest.mark.parametrize('wavelength, idl_offset',[
    (185.21*u.Angstrom, 16.992825*u.pixel),
    ([185.21]*u.Angstrom, 16.992825*u.pixel),
    (np.array([185.21])*u.Angstrom, 16.992825*u.pixel),
    (np.array([185.21, 195.119, 256.317, 258.375, 284.160])*u.Angstrom,
     [16.992825, 16.208033,  0.000239, -0.162755, -2.204928]*u.pixel)
])
def test_ccd_offset(wavelength, idl_offset):
    offset = ccd_offset(wavelength)
    idl_offset = np.atleast_1d(idl_offset)
    assert u.allclose(offset, idl_offset, atol=1e-5*u.pixel, rtol=0)

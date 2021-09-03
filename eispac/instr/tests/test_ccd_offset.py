import pytest

import numpy as np

from eispac.instr import ccd_offset


@pytest.mark.parametrize('wavelength, idl_offset',[
    (185.21, 16.992825),
    ([185.21], 16.992825),
    (np.array([185.21]), 16.992825),
    (np.array([185.21, 195.119, 256.317, 258.375, 284.160]),
     [16.992825, 16.208033,  0.000239, -0.162755, -2.204928])
])
def test_ccd_offset(wavelength, idl_offset):
    offset = ccd_offset(wavelength)
    idl_offset = np.atleast_1d(idl_offset)
    assert np.allclose(offset, idl_offset, atol=1e-5, rtol=0)

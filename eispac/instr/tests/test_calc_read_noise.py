import pytest
import numpy as np
from numpy.testing import assert_allclose

from eispac.instr.calc_read_noise import calc_read_noise

# Note: invalid inputs currently returns values of 0.0
# To-Do: refactor code to return "None" when given bad input values
def test_calc_read_noise_invalid_input():
    assert calc_read_noise(195.12) == 0.0
    assert calc_read_noise("195.12") == 0.0
    assert calc_read_noise({"wave": 195.12}) == 0.0
    assert calc_read_noise(None) == 0.0


@pytest.mark.parametrize("container_type", [list, tuple, np.array])
def test_calc_read_noise_container_types(container_type):
    wave_input = container_type([195.12, 256.32])
    rn = calc_read_noise(wave_input)

    assert isinstance(rn, np.ndarray)
    assert rn.shape == (2,)
    assert np.all(rn > 0)


def test_calc_read_noise_2d_and_3d_arrays():
    """Verify function preserves multidimensional array shapes (e.g. raster cubes)."""
    # 2D spectrum slice (e.g. n_slit_steps, n_wavelengths)
    wave_2d = np.ones((10, 20)) * 195.12
    rn_2d = calc_read_noise(wave_2d)
    assert rn_2d.shape == (10, 20)

    # 3D raster cube (e.g. n_pixels, n_steps, n_wavelengths)
    wave_3d = np.ones((5, 10, 20)) * 195.12
    rn_3d = calc_read_noise(wave_3d)
    assert rn_3d.shape == (5, 10, 20)


def test_calc_read_noise_formula():
    """
    Verify read noise calculation against the theoretical equation:
      read_noise_per_e = 14.427 electrons
      e_per_ph = (12398.5 / wavelength) / 3.65
      read_noise_counts = read_noise_per_e / e_per_ph
    """
    wave = np.array([180.0, 195.12, 270.0])
    expected_e_per_ph = (12398.5 / wave) / 3.65
    expected_rn = 14.427 / expected_e_per_ph

    rn = calc_read_noise(wave)
    assert isinstance(rn, np.ndarray)
    assert_allclose(rn, expected_rn, rtol=1e-7)


def test_calc_read_noise_proportionality():
    wave1 = np.array([100.0])
    wave2 = np.array([200.0])

    rn1 = calc_read_noise(wave1)
    rn2 = calc_read_noise(wave2)

    # Read noise in photon counts is directly proportional to wavelength
    # Doubling wavelength should double the read noise in photon counts
    assert_allclose(rn2, 2.0 * rn1, rtol=1e-7)




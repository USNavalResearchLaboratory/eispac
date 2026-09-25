import pytest
import numpy as np
import astropy.constants as const
from numpy.testing import assert_allclose

from eispac.instr.calc_velocity import calc_velocity

# Speed of light in km/s used by calc_velocity
C_KMS = const.c.to('km/s').value

# Note: invalid inputs currently returns values of 0.0
# To-Do: refactor code to return "None" when given bad input values
def test_calc_velocity_invalid_obs_wave():
    assert calc_velocity(195.119, 195.119) == 0.0
    assert calc_velocity("195.119", 195.119) == 0.0
    assert calc_velocity(None, 195.119) == 0.0


def test_calc_velocity_invalid_rest_wave():
    obs = [195.119, 195.120]
    assert calc_velocity(obs, [195.119]) == 0.0
    assert calc_velocity(obs, np.array(195.119)) == 0.0
    assert calc_velocity(obs, None) == 0.0
    # Malformed line ID strings
    assert calc_velocity(obs, "Fe_XII_195.119") == 0.0
    assert calc_velocity(obs, "Fe XII number") == 0.0


def test_calc_velocity_invalid_corr_method():
    obs = [195.119, 195.120]
    assert calc_velocity(obs, 195.119, corr_method='invalid_method') == 0.0
    assert calc_velocity(obs, 195.119, corr_method=123) == 0.0


def test_calc_velocity_zero_at_rest():
    vel = calc_velocity([195.119], 195.119, corr_method=None)
    assert_allclose(vel, [0.0], atol=1e-10)


def test_calc_velocity_basic_formula():
    rest_wave = 195.119
    obs_wave = np.array([195.119, 195.129, 195.109])
    expected_vel = C_KMS * (obs_wave - rest_wave) / rest_wave
    vel = calc_velocity(obs_wave, rest_wave, corr_method=None)
    assert isinstance(vel, np.ndarray)
    assert_allclose(vel, expected_vel, rtol=1e-7)


@pytest.mark.parametrize("input_type", [list, tuple, np.array])
def test_calc_velocity_input_container_types(input_type):
    obs_raw = [195.119, 195.125]
    obs_wave = input_type(obs_raw)
    vel = calc_velocity(obs_wave, 195.119, corr_method=None)

    assert isinstance(vel, np.ndarray)
    assert vel.shape == (2,)
    assert_allclose(vel[0], 0.0, atol=1e-10)


@pytest.mark.parametrize(
    "rest_wave_input, expected_rest",
    [
        (195.119, 195.119),
        (195, 195.0),
        ("195.119", 195.119),
        ("Fe XII 195.119", 195.119),
        ("Fe XXIV 192.03", 192.03),
    ],
)
def test_calc_velocity_rest_wave_formats(rest_wave_input, expected_rest):
    obs_wave = np.array([expected_rest, expected_rest + 0.01])
    expected_vel = C_KMS * (obs_wave - expected_rest) / expected_rest

    vel = calc_velocity(obs_wave, rest_wave_input, corr_method=None)
    assert_allclose(vel, expected_vel, rtol=1e-7)


def test_calc_velocity_corr_method_column():
    # 2D array of shape (n_pixels, n_steps) = (3, 2)
    obs_wave = np.array([
        [195.110, 195.120],
        [195.119, 195.130],
        [195.130, 195.140],
    ])
    rest_wave = 195.119

    # Uncorrected velocity
    uncorr_vel = C_KMS * (obs_wave - rest_wave) / rest_wave
    col_medians = np.median(uncorr_vel, axis=0)
    expected_vel = uncorr_vel - col_medians

    vel = calc_velocity(obs_wave, rest_wave, corr_method='column')
    assert vel.shape == (3, 2)
    assert_allclose(vel, expected_vel, rtol=1e-7)
    # Each column median should now be 0
    assert_allclose(np.median(vel, axis=0), [0.0, 0.0], atol=1e-10)


def test_calc_velocity_corr_method_image():
    obs_wave = np.array([
        [195.110, 195.120],
        [195.119, 195.130],
    ])
    rest_wave = 195.119

    uncorr_vel = C_KMS * (obs_wave - rest_wave) / rest_wave
    expected_vel = uncorr_vel - np.median(uncorr_vel)

    vel = calc_velocity(obs_wave, rest_wave, corr_method='image')
    assert vel.shape == (2, 2)
    assert_allclose(vel, expected_vel, rtol=1e-7)
    assert_allclose(np.median(vel), 0.0, atol=1e-10)



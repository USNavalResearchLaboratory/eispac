import pytest
import numpy as np
from numpy.testing import assert_allclose

from eispac.core.scale_guess import scale_guess

def test_scale_guess_wrong_num_params():
    x = np.linspace(195.0, 195.3, 31)
    y = np.ones(31)
    # Expected 3*1 + 0 = 3 params, but 4 provided
    param_guess = np.array([10.0, 195.15, 0.03, 0.0])

    with pytest.raises(SystemExit):
        scale_guess(x, y, param_guess, n_gauss=1, n_poly=0)
        scale_guess(x, y, param_guess, n_gauss=2, n_poly=1)


def test_scale_guess_single_gaussian_no_background():
    """
    Verify peak scaling for a single Gaussian with n_poly=0.
    Centroid and width should remain untouched, while peak is set to y at centroid index.
    """
    x = np.linspace(195.0, 195.3, 31)
    # Synthetic Gaussian: peak = 150.0, center = 195.15, width = 0.03
    y = 150.0 * np.exp(-((x - 195.15) ** 2) / (2 * 0.03 ** 2))

    # Initial guess with incorrect peak
    param_guess = np.array([10.0, 195.15, 0.03])

    scaled_param = scale_guess(x, y, param_guess, n_gauss=1, n_poly=0)

    assert len(scaled_param) == 3 # num of params unchanged
    # Peak should be scaled to data maximum (~150)
    assert_allclose(scaled_param[0], 150.0, rtol=1e-3)
    # Centroid and width should remain unchanged
    assert_allclose(scaled_param[1], 195.15)
    assert_allclose(scaled_param[2], 0.03)


def test_scale_guess_multigaussian():
    """
    Verify scaling across multiple Gaussian components simultaneously.
    """
    bg_level = 10.0
    x = np.linspace(195.0, 195.5, 51)
    # Two peaks: line 1 at 195.15 (peak=100), line 2 at 195.35 (peak=50)
    y = (bg_level
        + 100.0 * np.exp(-((x - 195.15) ** 2) / (2 * 0.02 ** 2))
        + 50.0 * np.exp(-((x - 195.35) ** 2) / (2 * 0.02 ** 2))
    )

    # Initial guess layout: [p1, c1, w1, p2, c2, w2, bg]
    param_guess = np.array([1.0, 195.15, 0.02, 1.0, 195.35, 0.02, 2.0])
    scaled_param = scale_guess(x, y, param_guess, n_gauss=2, n_poly=1)

    assert len(scaled_param) == 7
    # Background parameter (index 6) should be scaled to data background (~10.0)
    assert_allclose(scaled_param[6], bg_level, rtol=1e-2)
    # Peak 1 (index 0) and Peak 2 (index 3) scaled after subtracting background
    assert_allclose(scaled_param[0], 100.0, rtol=1e-2)
    assert_allclose(scaled_param[3], 50.0, rtol=1e-2)
    # Centroids and widths unchanged
    assert_allclose(scaled_param[[1, 2, 4, 5]], [195.15, 0.02, 195.35, 0.02])


def test_scale_guess_negative_peak_clipping():
    """
    When y[centroid] - background is negative, the peak should be clipped to 0.0.
    """
    x = np.linspace(195.0, 195.3, 31)
    y = np.full_like(x, 20.0)
    # Dip at "centroid" location below background level
    y[15] = 5.0

    # Guess expecting a peak at index 15
    param_guess = np.array([10.0, x[15], 0.03, 20.0])
    scaled_param = scale_guess(x, y, param_guess, n_gauss=1, n_poly=1)

    assert scaled_param[0] == 0.0
import numpy as np
import pytest
from numpy.testing import assert_allclose

import eispac


def test_multigaussian_invalid_param_size():
    with pytest.raises(SystemExit):
        func_vals = eispac.multigaussian([1], [1,2,3], n_gauss=1, n_poly=0)


def test_deviates_invalid_param_size_debug():
    with pytest.raises(SystemExit):
        deviates = eispac.multigaussian_deviates([1], x=[1,2,3], y=[1,4,9],
                                                 n_gauss=1, n_poly=0, debug=True)


def test_deviates_invalid_param_size():
    deviates = eispac.multigaussian_deviates([1], x=[1,2,3], y=[1,4,9],
                                             n_gauss=1, n_poly=0)
    assert deviates[0] == -3

# ----------------------------------------------------------------------
# 1. multigaussian() Mathematical Output Tests
# ----------------------------------------------------------------------
#TO-DO: Figure out what the point of some of these tests are!
#TO-DO: Update with proper and consistent x resolution
def test_multigaussian_single_comp():
    """Verify single Gaussian evaluation at peak, 1-sigma, and symmetric points."""
    # x = np.linspace(195.0, 195.3, 31)
    # param = [100.0, 195.119, 0.03] #[peak, centroid, width]
    x = np.array([195.0, 195.119, 195.238]) # WHY ONLY 3 DATA POINT?!?
    params = [100.0, 195.119, 0.05] #[peak, centroid, width]

    y_eval = eispac.multigaussian(params, x, n_gauss=1, n_poly=0)

    # At centroid (x[1]), f(x) must equal peak
    assert_allclose(y_eval[1], params[0], rtol=1e-7)

    # WHY ARE WE COMPUTING ADDITIONAL GAUSSIAN PROFILES?!?!
    # Symmetric points equidistant from centroid must yield identical values
    dx = 0.03
    y_left = eispac.multigaussian(params, np.array([params[1] - dx]), n_gauss=1, n_poly=0)
    y_right = eispac.multigaussian(params, np.array([params[1] + dx]), n_gauss=1, n_poly=0)
    assert_allclose(y_left, y_right, rtol=1e-7)

    # Verify 1-sigma value: f(centroid +/- width) == peak * exp(-0.5)
    y_1sigma = eispac.multigaussian(params, np.array([params[1] + params[2]]), n_gauss=1, n_poly=0)
    assert_allclose(y_1sigma, params[0] * np.exp(-0.5), rtol=1e-7)


def test_multigaussian_multi_comp_with_background():
    """Verify multi-component Gaussian combined with polynomial background."""
    x = np.linspace(194.0, 196.0, 50)
    # Params: [p1, c1, w1, p2, c2, w2, c0, c1]
    # Poly terms are evaluated via np.polynomial.polyval (lowest order first: c0 + c1*x)
    params = [80.0, 194.8, 0.04, 40.0, 195.3, 0.03, 10.0, 0.5]
    n_gauss = 2
    n_poly = 2

    y_eval = eispac.multigaussian(params, x, n_gauss=n_gauss, n_poly=n_poly)

    # Check order of background terms
    assert y_eval[0] == params[-2] + params[-1]*x[0]

    # Expected value: Gaussian 1 + Gaussian 2 + Polynomial
    g1 = 80.0 * np.exp(-0.5 * ((x - 194.8) / 0.04) ** 2)
    g2 = 40.0 * np.exp(-0.5 * ((x - 195.3) / 0.03) ** 2)
    poly = 10.0 + 0.5 * x
    expected = g1 + g2 + poly

    assert_allclose(y_eval, expected, rtol=1e-10)


# ----------------------------------------------------------------------
# 2. multigaussian_deviates() Residuals & Data Masking Tests
# ----------------------------------------------------------------------
def test_multigaussian_deviates_valid_data():
    """Verify standard residuals calculation (y - model) / error."""
    x = np.linspace(194.5, 195.5, 20)
    param = [100.0, 195.0, 0.05, 5.0]
    model = eispac.multigaussian(param, x, n_gauss=1, n_poly=1)

    # Generate synthetic observed y and error
    err = np.full_like(x, 2.0)
    y_obs = model + 1.5
    expected_deviates = (y_obs - model) / err

    status, deviates = eispac.multigaussian_deviates(
        param, x=x, y=y_obs, error=err, n_gauss=1, n_poly=1
    )

    assert status == 0
    assert_allclose(deviates, expected_deviates, rtol=1e-10)


def test_multigaussian_deviates_masked_data():
    """Verify points with negative error are zeroed out in output deviates."""
    x = np.linspace(194.5, 195.5, 10)
    param = [100.0, 195.0, 0.05]
    y_obs = np.ones(10) * 50.0
    err = np.array([1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    status, deviates = eispac.multigaussian_deviates(
        param, x=x, y=y_obs, error=err.copy(), n_gauss=1, n_poly=0
    )

    assert status == 0
    # Indices 2 and 4 had error < 0, so their deviates must be 0.0
    assert deviates[2] == 0.0
    assert deviates[4] == 0.0
    # Good points should not be zero
    assert deviates[0] != 0.0


def test_multigaussian_deviates_all_bad_data():
    """Verify return status -2 when all data points have negative errors."""
    x = np.array([195.0, 195.1, 195.2])
    y_obs = np.array([10.0, 20.0, 10.0])
    err = np.array([-1.0, -1.0, -1.0])
    param = [100.0, 195.1, 0.05]

    # Non-debug mode returns status -2
    status, deviates = eispac.multigaussian_deviates(
        param, x=x, y=y_obs, error=err.copy(), n_gauss=1, n_poly=0, debug=False
    )
    assert status == -2
    assert_allclose(deviates, np.full(3, -1.0))

    # Debug mode raises SystemExit
    with pytest.raises(SystemExit):
        eispac.multigaussian_deviates(
            param, x=x, y=y_obs, error=err.copy(), n_gauss=1, n_poly=0, debug=True
        )
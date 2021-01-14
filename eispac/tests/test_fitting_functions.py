import numpy as np
import pytest
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

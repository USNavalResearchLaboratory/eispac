import pytest
import numpy as np

import eispac


@pytest.fixture
def empty_template():
    return eispac.EISFitTemplate()


def test_empty_is_fit_template(empty_template):
    assert isinstance(empty_template, eispac.EISFitTemplate)


def empty_template_has_repr(empty_template):
    assert isinstance(empty_template.__repr__(), str)

def empty_template_has_template(empty_template):
    assert isinstance(empty_template.template, dict)
    assert all(k in empty_template.template for k in ['n_gauss', 'n_poly', 'fit', 'line_ids', 'wmin', 'wmax'])

def empty_template_has_parinfo(empty_template):
    assert isinstance(empty_template.parinfo, list)
    for f in empty_template.parinfo:
        assert all(k in f for k in ['fixed', 'limited', 'limits', 'tied', 'value'])

def empty_template_has_funcinfo(empty_template):
    assert isinstance(empty_template.funcinfo, list)
    for f in empty_template.funcinfo:
        assert all(k in f for k in ['func', 'name', 'n_params'])

def empty_template_parameters_length(empty_template):
    n_params = sum([f['n_params'] for f in empty_template.funcinfo])
    assert n_params == len(empty_template.parinfo)


def empty_template_central_wave(empty_template):
    assert isinstance(empty_template.central_wave, float)

import pytest

import eispac
from eispac.templates import fit_template_filenames


@pytest.fixture(params=fit_template_filenames())
def test_template(request):
    return eispac.EISFitTemplate.read_template(request.param)


def test_template_is_template(test_template):
    assert isinstance(test_template, eispac.EISFitTemplate)


def test_template_has_repr(test_template):
    assert isinstance(test_template.__repr__(), str)


def test_template_has_funcinfo(test_template):
    assert isinstance(test_template.funcinfo, list)
    for f in test_template.funcinfo:
        assert all(k in f for k in ['func', 'name', 'n_params'])


def test_template_has_parinfo(test_template):
    assert isinstance(test_template.parinfo, list)
    for f in test_template.parinfo:
        assert all(k in f for k in ['fixed', 'limited', 'limits', 'tied', 'value'])


def test_template_parameters_length(test_template):
    n_params = sum([f['n_params'] for f in test_template.funcinfo])
    assert n_params == len(test_template.parinfo)


def test_template_central_wave(test_template):
    assert isinstance(test_template.central_wave, float)

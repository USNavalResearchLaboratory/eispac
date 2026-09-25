import pathlib
import warnings

import pytest

import eispac
from eispac.data import fit_template_filenames


# Note: Since the addition of custom templates in EISPAC v0.95 (2024-03-05),
#       eispac now validates the contents of a template after loading it.
#       Therefore, we now only need to test that a given template can
#       be read without warnings being printed to stderr.
#       Tests for the actual validation routines are in "test_eisfittemplate" 
@pytest.mark.parametrize('filename', [42, '/not/a/valid/path', pathlib.Path('/not/a/valid/path')])
def test_bad_filename_returns_none(filename):
    with pytest.warns(UserWarning, match='Error: Template filepath'):
        assert eispac.EISFitTemplate.read_template(filename) is None


def test_read_template_hdf5(test_template_filepath):
    hdf5_template = eispac.read_template(test_template_filepath)
    assert isinstance(hdf5_template, eispac.EISFitTemplate)


def test_read_template_toml(test_toml_template_filepath):
    toml_template = eispac.read_template(test_toml_template_filepath)
    assert isinstance(toml_template, eispac.EISFitTemplate)


@pytest.mark.parametrize('template_filepath', fit_template_filenames())
def test_read_template_premade(capsys, template_filepath):
    test_template = eispac.read_template(template_filepath)

    # Check that no warnings were written to stderr
    captured = capsys.readouterr()
    assert len(captured.err) == 0

import numpy as np
import numpy.testing as npt
import pytest

from eispac.core.eisfittemplate import EISFitTemplate


@pytest.fixture
def empty_template():
    return EISFitTemplate()


def test_empty_is_fit_template(empty_template):
    assert isinstance(empty_template, EISFitTemplate)
    assert isinstance(empty_template.__repr__(), str)


def test_empty_template_has_template(empty_template):
    assert isinstance(empty_template.template, dict)
    assert all(k in empty_template.template 
               for k in ['n_gauss', 'n_poly', 'line_ids', 'wmin', 'wmax', 'fit'])


def test_empty_template_has_parinfo(empty_template):
    assert isinstance(empty_template.parinfo, list)
    for f in empty_template.parinfo:
        assert all(k in f for k in ['value', 'fixed', 'limited', 'limits', 'tied'])


def test_empty_template_has_funcinfo(empty_template):
    assert isinstance(empty_template.funcinfo, list)
    for f in empty_template.funcinfo:
        assert all(k in f for k in ['func', 'name', 'n_params'])


def test_empty_template_parameters_length(empty_template):
    n_params = sum([f['n_params'] for f in empty_template.funcinfo])
    assert n_params == len(empty_template.parinfo)


def empty_template_central_wave(empty_template):
    assert isinstance(empty_template.central_wave, float)


# ======================================================================
# Error and Exception Tests for _validate_parinfo_list
# ======================================================================
def test_parinfo_invalid_datatype(capsys):
    """Passing a non-list/dict parinfo logs an error to stderr and initializes defaults."""
    tmpl = EISFitTemplate(parinfo=42)
    captured = capsys.readouterr()

    assert "Error: Invalid datatype for \"parinfo\"" in captured.err
    assert isinstance(tmpl.parinfo, list)
    assert len(tmpl.parinfo) > 0


def test_parinfo_missing_value_key(capsys):
    """Parinfo dict without a 'value' key logs a warning and estimates parameter count."""
    dict_parinfo = {"fixed": [0, 0, 1],
                    "step": [0.1, 0.1, 0.1]}
    
    tmpl = EISFitTemplate(template={"n_gauss": 1, "n_poly": 0}, parinfo=dict_parinfo)
    captured = capsys.readouterr()

    assert "Warning: No initial values given in parinfo!" in captured.err
    assert len(tmpl.parinfo) == 3
    # Missing 'value' should be populated with default 1.0
    for p in tmpl.parinfo:
        assert p["value"] == 1.0


def test_parinfo_unsupported_key_ignored(capsys):
    """Unsupported parinfo keys log a warning and are excluded."""
    dict_parinfo = {"value": [10.0, 195.12, 0.04],
                    "invalid_custom_key": [1, 2, 3]}
    
    tmpl = EISFitTemplate(template={"n_gauss": 1, "n_poly": 0}, parinfo=dict_parinfo)
    captured = capsys.readouterr()

    assert "is not a valid parinfo key and will be ignored" in captured.err
    for p in tmpl.parinfo:
        assert "invalid_custom_key" not in p


def test_parinfo_type_conversion_failure(capsys):
    """Values that cannot be cast to expected types log a warning and skip the key."""
    dict_parinfo = {"value": [10.0, 195.12, 0.04],
                    # 'fixed' expects an int, pass unparseable strings
                    "fixed": ["cannot_convert", 0, 1]}
    
    tmpl = EISFitTemplate(template={"n_gauss": 1, "n_poly": 0}, parinfo=dict_parinfo)
    captured = capsys.readouterr()

    assert "There was an issue loading the fixed key for parameter index 0" in captured.err
    # Parameter 0 should fallback to default fixed=0
    assert tmpl.parinfo[0]["fixed"] == 0
    assert tmpl.parinfo[2]["fixed"] == 1


@pytest.mark.parametrize("invalid_shape", [[1], [1, 0, 0], [],])
def test_parinfo_invalid_limits_and_limited_length(capsys, invalid_shape):
    """'limited' or 'limits' with length != 2 logs a warning and resets to default."""
    parinfo = [{"value": 10.0, "limited": invalid_shape},
               {"value": 195.12, "limits": invalid_shape},
               {"value": 0.04}]
    
    tmpl = EISFitTemplate(template={"n_gauss": 1, "n_poly": 0}, parinfo=parinfo)
    captured = capsys.readouterr()

    assert 'Warning: incorrect length of "limited" and/or "limits"' in captured.err
    # Validated structures must be replaced by standard 2-element arrays
    assert len(tmpl.parinfo[0]["limited"]) == 2
    assert len(tmpl.parinfo[1]["limits"]) == 2
    npt.assert_array_equal(tmpl.parinfo[0]["limited"], np.array([1, 0], dtype="int16"))
    npt.assert_array_equal(tmpl.parinfo[1]["limits"], np.zeros(2))

# ----------------------------------------------------------------------
# Tests for _validate_parinfo_list
# ----------------------------------------------------------------------
def test_parinfo_basic_input():
    """Verify parinfo is read and missing keys are added"""
    # 1 Gaussian (3 params) + 0 poly = 3 params
    list_parinfo = [{"value": 100.0}, {"VALUE": 195.12}, {"VaLuE": 0.05}]
    
    fit_tmpl = EISFitTemplate(template={"n_gauss": 1, "n_poly": 0}, parinfo=list_parinfo)

    assert len(fit_tmpl.parinfo) == 3
    # Note: these test also confirm the input keys are case insensitive
    for p in fit_tmpl.parinfo:
        assert "fixed" in p and p["fixed"] == 0
        assert "limited" in p and isinstance(p["limited"], np.ndarray)
        assert "limits" in p and isinstance(p["limits"], np.ndarray)
        assert "tied" in p
        assert "value" in p
        assert "VALUE" not in p
        assert "VaLuE" not in p

    # Check default arrays
    npt.assert_array_equal(fit_tmpl.parinfo[0]["limited"], 
                           np.array([1, 0], dtype="int16"))
    npt.assert_array_equal(fit_tmpl.parinfo[0]["limits"], np.zeros(2))


def test_parinfo_dict_of_lists_conversion():
    """Verify conversion from a dict-of-lists/arrays to a list-of-dicts."""
    dict_parinfo = {"value": [10.0, 195.12, 0.04, 1.5],
                    "fixed": [0, 0, 0, 1],
                    "limited": [[1, 0], [1, 1], [1, 0], [0, 0]],
                    "limits": [[0.0, 0.0], [195.0, 195.3], [0.0, 0.0], [0.0, 0.0]]}
    
    fit_tmpl = EISFitTemplate(template={"n_gauss": 1, "n_poly": 1}, parinfo=dict_parinfo)

    assert isinstance(fit_tmpl.parinfo, list)
    assert len(fit_tmpl.parinfo) == 4
    assert fit_tmpl.parinfo[0]["value"] == 10.0
    assert fit_tmpl.parinfo[1]["value"] == 195.12
    assert fit_tmpl.parinfo[3]["fixed"] == 1
    npt.assert_array_equal(fit_tmpl.parinfo[1]["limits"], np.array([195.0, 195.3]))


def test_parinfo_custom_mpfit_keys_preserved():
    """Verify that additional MPFIT keys like 'parname', 'step', 'tied' are preserved."""
    custom_parinfo = [{"value": 10.0, "parname": "Peak", "step": 0.1},
                      {"value": 195.12, "parname": "Centroid", "tied": "p[4] - 0.5"},
                      {"value": 0.04, "parname": "Width", "mpmaxstep": 0.01}]
    
    fit_tmpl = EISFitTemplate(template={"n_gauss": 1, "n_poly": 0}, 
                              parinfo=custom_parinfo)

    assert fit_tmpl.parinfo[0]["parname"] == "Peak"
    assert fit_tmpl.parinfo[0]["step"] == 0.1
    assert "p[4] - 0.5" in str(fit_tmpl.parinfo[1]["tied"])
    assert fit_tmpl.parinfo[2]["mpmaxstep"] == 0.01

# ======================================================================
# Error and Exception Tests for _validate_template_dict
# ======================================================================
def test_template_invalid_datatype(capsys):
    """Passing a non-dict template logs an error and uses an empty dict."""
    tmpl = EISFitTemplate(template="invalid_template_type")
    captured = capsys.readouterr()

    assert "Error: Invalid datatype for \"template\"" in captured.err
    assert isinstance(tmpl.template, dict)


@pytest.mark.parametrize("invalid_ngauss", [0, -1, -5, "one", None, []])
def test_template_invalid_n_gauss(capsys, invalid_ngauss):
    """n_gauss <= 0 or invalid datatype logs an error and sets n_gauss to -1."""
    tmpl = EISFitTemplate()
    tmpl.template["n_gauss"] = invalid_ngauss
    tmpl._validate_template_dict()
    captured = capsys.readouterr()

    assert "Error: Invalid value for n_gauss" in captured.err


@pytest.mark.parametrize("invalid_npoly", [-1, -3, "zero", None, {}])
def test_template_invalid_n_poly(capsys, invalid_npoly):
    """n_poly < 0 or invalid datatype logs an error and sets n_poly to -1."""
    tmpl = EISFitTemplate()
    tmpl.template["n_poly"] = invalid_npoly
    tmpl._validate_template_dict()
    captured = capsys.readouterr()

    assert "Error: Invalid value for n_poly" in captured.err


def test_template_invalid_wmin_and_wmax(capsys):
    """Non-numeric wmin and wmax log errors and reset to defaults (170.0, 292.0)."""
    tmpl = EISFitTemplate(template={"n_gauss": 1, "n_poly": 0, 
                                    "wmin": "not_a_float", "wmax": [292.0]},
                          parinfo=[{"value": 10.0}, {"value": 195.12}, {"value": 0.04}])
    captured = capsys.readouterr()

    assert "Error: Invalid datatype for wmin" in captured.err
    assert "Error: Invalid datatype for wmax" in captured.err
    assert tmpl.template["wmin"] == 170.0
    assert tmpl.template["wmax"] == 292.0


def test_template_invalid_component(capsys):
    """Non-numeric component logs an error and defaults to 1."""
    tmpl = EISFitTemplate(template={"n_gauss": 1, "n_poly": 0, "component": "comp_1"},
                          parinfo=[{"value": 10.0}, {"value": 195.12}, {"value": 0.04}])
    captured = capsys.readouterr()

    assert "Error: Invalid datatype for component" in captured.err
    assert tmpl.template["component"] == 1


def test_template_parameter_count_mismatch(capsys):
    """When n_gauss and n_poly do not match parinfo length, log warning and auto-adjust."""
    # 4 parameters given, but template specifies n_gauss=2, n_poly=1 (requires 7)
    # parinfo is assumed to be correct
    parinfo = [{"value": 10.0}, {"value": 195.12}, {"value": 0.04}, {"value": 1.0}]
    tmpl = EISFitTemplate(template={"n_gauss": 2, "n_poly": 1}, parinfo=parinfo)
    captured = capsys.readouterr()

    assert "Warning: the values of n_gauss and n_poly do not match" in captured.err
    # 4 parameters total -> 1 Gaussian (3 params) + 1 Poly (1 param)
    assert tmpl.template["n_gauss"] == 1
    assert tmpl.template["n_poly"] == 1
    assert len(tmpl.template["line_ids"]) == 1


def test_template_default_line_ids():
    """When fewer line_ids are provided than n_gauss, missing IDs are filled with defaults."""
    parinfo = [{"value": 10.0}, {"value": 195.12}, {"value": 0.04},  # Gauss 1
               {"value": 5.0},  {"value": 195.40}, {"value": 0.03}]  # Gauss 2

    # n_gauss = 2, but only 1 line_id provided
    tmpl = EISFitTemplate(template={"n_gauss": 2, "n_poly": 0, 
                                    "line_ids": ["Fe XII 195.119"]},
                          parinfo=parinfo)

    assert len(tmpl.template["line_ids"]) == 2
    assert tmpl.template["line_ids"][0] == "Fe XII 195.119"
    # Second item should be auto-populated with placeholder
    assert "unknown" in tmpl.template["line_ids"][1].lower()
    assert "195.400" in tmpl.template["line_ids"][1].lower()

# ----------------------------------------------------------------------
# Tests for _validate_template_dict
# ----------------------------------------------------------------------
def test_template_basic_input():
    """Verify template dict is loaded and missing keys are added"""

    fit_tmpl = EISFitTemplate(template={"n_gauss": 1, "n_poly": 0, 
                                        "WMIN": 194.5, "wmax": 195.5},
                              parinfo=[{"value": 10.0}, {"value": 195.12}, 
                                       {"value": 0.04}])

    # Check that wmin and wmax are unchanged
    assert "WMIN" not in fit_tmpl.template.keys() # Check that keys are lowercased
    assert fit_tmpl.template["wmin"] == 194.5
    assert fit_tmpl.template["wmax"] == 195.5

    # Check that num of params matches parinfo
    assert fit_tmpl.template["n_gauss"] == 1
    assert fit_tmpl.template["n_poly"] == 0
    assert len(fit_tmpl.parinfo) == 3

    # Check that missing template keys are added
    assert 'fit' in fit_tmpl.template.keys()
    assert 'line_ids' in fit_tmpl.template.keys()
    npt.assert_array_equal(fit_tmpl.template["fit"], np.array([10.0, 195.12, 0.04]))


# ----------------------------------------------------------------------
# Tests for _check_param_values
# ----------------------------------------------------------------------
def test_check_param_values_wmin_greater_than_wmax(capsys):
    """Verify warning when wmin is greater than wmax."""
    parinfo = [{"value": 10.0}, {"value": 195.0}, {"value": 0.04}]
    template = {"n_gauss": 1, "n_poly": 0, "wmin": 196.0, "wmax": 194.0}

    tmpl = EISFitTemplate(template=template, parinfo=parinfo)
    captured = capsys.readouterr()

    assert "Warning: wmin" in captured.err


def test_check_param_values_large_wavelength_range(capsys):
    """Verify warning when wmax - wmin exceeds 10 Angstroms."""
    parinfo = [{"value": 10.0}, {"value": 190.0}, {"value": 0.04}]
    template = {"n_gauss": 1, "n_poly": 0, "wmin": 180.0, "wmax": 195.0}

    tmpl = EISFitTemplate(template=template, parinfo=parinfo)
    captured = capsys.readouterr()

    assert "Warning: wmin and wmax define" in captured.err


def test_check_param_values_centroid_outside_range(capsys):
    """Verify warning when a Gaussian centroid is outside [wmin, wmax]."""
    parinfo = [{"value": 10.0}, {"value": 190.0}, {"value": 0.04},  # Centroid 0
               {"value": 15.0}, {"value": 195.0}, {"value": 0.04}]  # Centroid 1
    template = {"n_gauss": 2, "n_poly": 0, "wmin": 194.0, "wmax": 196.0}

    tmpl = EISFitTemplate(template=template, parinfo=parinfo)
    captured = capsys.readouterr()

    assert "Warning: centroid for Gaussian 0" in captured.err
    assert "outside of the fitting wavelength range" in captured.err
    assert "Warning: centroid for Gaussian 1" not in captured.err


def test_check_param_values_valid_inputs_no_warnings(capsys):
    """Verify no warnings are emitted when bounds and centroids are valid."""
    parinfo = [{"value": 10.0}, {"value": 195.12}, {"value": 0.04}]
    template = {"n_gauss": 1, "n_poly": 0, "wmin": 194.5, "wmax": 195.5}

    tmpl = EISFitTemplate(template=template, parinfo=parinfo)
    captured = capsys.readouterr()

    assert "Warning" not in captured.err
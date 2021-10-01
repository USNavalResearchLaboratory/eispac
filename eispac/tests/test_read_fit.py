import eispac


def test_invalid_filepath_type():
    fit_res = eispac.read_fit(42)
    assert fit_res is None

def test_missing_file():
    fit_res = eispac.read_fit('non-existent_file.fit.h5')
    assert fit_res is None

def test_read_fit_string_filepath(test_fit_filepath):
    fit_res = eispac.read_fit(test_fit_filepath)
    assert isinstance(fit_res, eispac.EISFitResult)

def test_read_fit_pathlib_filepath(test_fit_filepath):
    import pathlib
    path_obj = pathlib.Path(test_fit_filepath)
    fit_res = eispac.read_fit(path_obj)
    assert isinstance(fit_res, eispac.EISFitResult)

def test_read_fit_verbose(test_fit_filepath):
    fit_res = eispac.read_fit(test_fit_filepath, verbose=True)
    assert isinstance(fit_res, eispac.EISFitResult)

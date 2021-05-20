import pytest
import eispac

# test_tmplt_filepath = '../templates/eis_template_dir/fe_12_192_394.1c.template.h5'
test_tmplt_filepath = eispac.templates.template_dir+'/fe_12_192_394.1c.template.h5'

def test_invalid_filepath_type():
    tmplt = eispac.read_template(42)
    assert tmplt is None

def test_missing_file():
    tmplt = eispac.read_template('non-existent_file.head.h5')
    assert tmplt is None

def test_read_template_str_filepath():
    tmplt = eispac.read_template(test_tmplt_filepath, quiet=True)
    assert isinstance(tmplt, eispac.EISFitTemplate)

def test_read_template_pathlib_filepath():
    import pathlib
    path_obj = pathlib.Path(test_tmplt_filepath)
    tmplt = eispac.read_template(path_obj, quiet=True)
    assert isinstance(tmplt, eispac.EISFitTemplate)

def test_print_parinfo():
    tmplt = eispac.read_template(test_tmplt_filepath, quiet=True)
    tmplt.print_parinfo()
    assert isinstance(tmplt, eispac.EISFitTemplate)

def test_create_funcinfo():
    tmplt = eispac.read_template(test_tmplt_filepath)
    assert len(tmplt.funcinfo) == 2

import pytest

from eispac.data import get_test_filepath
from eispac.templates import get_template_filepath


@pytest.fixture
def test_data_filepath():
    return get_test_filepath('eis_20210306_064444.data.h5')


@pytest.fixture
def test_head_filepath():
    return get_test_filepath('eis_20210306_064444.head.h5')


@pytest.fixture
def test_fit_filepath():
    return get_test_filepath('eis_20210306_064444.fe_12_192_394.1c-0.fit.h5')


@pytest.fixture
def test_template_filepath():
    return get_template_filepath('fe_12_192_394.1c.template.h5')
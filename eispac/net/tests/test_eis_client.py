import pytest
from sunpy.net import Fido, attrs as a
import eispac.net


@pytest.fixture
def eis_query():
    time = a.Time('2022-03-29 22:21:00','2022-03-29 23:21:00')
    instr = a.Instrument('EIS')
    obs = a.Physobs.intensity
    source = a.Source('Hinode')
    provider = a.Provider('NRL')
    level = a.Level('1')
    return time, instr, obs, source, provider, level


@pytest.mark.remote_data
def test_search_all_types(eis_query):
    q = Fido.search(*eis_query)
    assert len(q) == 1
    assert len(q[0]) == 3
    assert q[0,0]['url'] == 'https://eis.nrl.navy.mil/level1/hdf5/2022/03/29/eis_20220329_222113.data.h5'


@pytest.mark.remote_data
def test_search_fits_only(eis_query):
    q = Fido.search(*eis_query, a.eispac.FileType('FITS'))
    assert len(q) == 1
    assert len(q[0]) == 1
    assert q[0,0]['url'] == 'https://eis.nrl.navy.mil/level1/fits/2022/03/29/eis_er_20220329_222113.fits'


@pytest.mark.parametrize('file_type, file_url', [
    ('FITS', 'https://eis.nrl.navy.mil/level1/fits/2022/03/29/eis_er_20220329_222113.fits'),
    ('HDF5 data', 'https://eis.nrl.navy.mil/level1/hdf5/2022/03/29/eis_20220329_222113.data.h5'),
    ('HDF5 header', 'https://eis.nrl.navy.mil/level1/hdf5/2022/03/29/eis_20220329_222113.head.h5')
])
@pytest.mark.remote_data
def test_search_individual_filetypes(eis_query, file_type, file_url):
    q = Fido.search(*eis_query, a.eispac.FileType(file_type))
    assert len(q) == 1
    assert len(q[0]) == 1
    assert q[0,0]['url'] == file_url
    assert q[0,0]['FileType'] == file_type


@pytest.mark.remote_data
def test_combined_hdf5_search(eis_query):
    q = Fido.search(*eis_query,
                    a.eispac.FileType('HDF5 data') | a.eispac.FileType('HDF5 header'))
    assert len(q) == 1
    assert len(q[0]) == 2
    assert q[0,0]['FileType'] == 'HDF5 data'
    assert q[0,1]['FileType'] == 'HDF5 header'


def test_registered_attrs():
    attr_names = ['fits', 'data_h5', 'head_h5']
    for an in attr_names:
        assert hasattr(a.eispac.FileType, an)

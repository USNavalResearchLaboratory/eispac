import urllib
import pytest
from sunpy.net import Fido, attrs as a
import eispac.net

# Check status of remote servers so we know when to skip the tests
try:
    with urllib.request.urlopen('https://eis.nrl.navy.mil/level1/') as nrl_site:
        nrl_status = nrl_site.status
except:
    nrl_status = -1
nrl_unavailable = nrl_status < 200 or nrl_status >= 300

try:
    with urllib.request.urlopen('https://vsolar.mssl.ucl.ac.uk/eispac/hdf5/') as mssl_site:
        mssl_status = nrl_site.status
except:
    mssl_status = -1
mssl_unavailable = mssl_status < 200 or mssl_status >= 300

@pytest.fixture
def eis_query():
    time = a.Time('2022-03-29 22:21:00','2022-03-29 23:21:00')
    instr = a.Instrument('EIS')
    obs = a.Physobs.intensity
    source = a.Source('Hinode')
    provider = a.Provider('NRL')
    level = a.Level('1')
    return time, instr, obs, source, provider, level

def test_registered_attrs():
    attr_names = ['data_h5', 'head_h5', 'fits', 'any']
    for an in attr_names:
        assert hasattr(a.eispac.FileType, an)

@pytest.mark.skipif(nrl_unavailable, reason="NRL server issues")
def test_search_default(eis_query):
    q = Fido.search(*eis_query)
    assert len(q) == 1
    assert len(q[0]) == 1
    assert q[0,0]['url'] == 'https://eis.nrl.navy.mil/level1/hdf5/2022/03/29/eis_20220329_222113.data.h5'

@pytest.mark.skipif(nrl_unavailable, reason="NRL server issues")
def test_search_all_types(eis_query):
    q = Fido.search(*eis_query, a.eispac.FileType('Any'))
    assert len(q) == 1
    assert len(q[0]) == 3
    assert q[0,0]['url'] == 'https://eis.nrl.navy.mil/level1/hdf5/2022/03/29/eis_20220329_222113.data.h5'

@pytest.mark.parametrize('file_type, file_url', [
    ('FITS', 'https://eis.nrl.navy.mil/level1/fits/2022/03/29/eis_l1_20220329_222113.fits'),
    ('HDF5 data', 'https://eis.nrl.navy.mil/level1/hdf5/2022/03/29/eis_20220329_222113.data.h5'),
    ('HDF5 header', 'https://eis.nrl.navy.mil/level1/hdf5/2022/03/29/eis_20220329_222113.head.h5')
])

@pytest.mark.skipif(nrl_unavailable, reason="NRL server issues")
def test_search_individual_filetypes(eis_query, file_type, file_url):
    q = Fido.search(*eis_query, a.eispac.FileType(file_type))
    assert len(q) == 1
    assert len(q[0]) == 1
    assert q[0,0]['url'] == file_url
    assert q[0,0]['FileType'] == file_type

@pytest.mark.skipif(nrl_unavailable, reason="NRL server issues")
def test_combined_hdf5_search(eis_query):
    q = Fido.search(*eis_query,
                    a.eispac.FileType('HDF5 data') | a.eispac.FileType('HDF5 header'))
    assert len(q) == 2
    assert len(q[0]) == 1
    assert len(q[1]) == 1
    assert q[0,0]['FileType'] == 'HDF5 data'
    assert q[1,0]['FileType'] == 'HDF5 header'

@pytest.mark.skipif(mssl_unavailable, reason="MSSL server issues")
def test_search_mssl(eis_query):
    q = Fido.search(*eis_query[:4], a.Provider('MSSL'))
    assert len(q) == 1
    assert len(q[0]) == 1
    assert q[0,0]['url'] == 'https://vsolar.mssl.ucl.ac.uk/eispac/hdf5/2022/03/29/eis_20220329_222113.data.h5'

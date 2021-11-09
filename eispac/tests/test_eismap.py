"""
Tests for `eispac.EISMap`
"""
import copy
from textwrap import dedent

import numpy as np
import pytest
from astropy.io import fits
import astropy.time
import astropy.units as u
from astropy.visualization import ImageNormalize, AsinhStretch
import sunpy
import sunpy.map
from sunpy.time import parse_time

# This registers the EISMap map class
import eispac  # NOQA


@pytest.fixture
def test_header():
    raw_header = dedent("""\
        SIMPLE  =                    T / conforms to FITS standard
        BITPIX  =                  -64 / array data type
        NAXIS   =                    2 / number of array dimensions
        NAXIS1  =                   60
        NAXIS2  =                  512
        EXTEND  =                    T
        DATE_OBS= '2010-07-23T14:32:10.000'
        DATE_BEG= '2010-07-23T14:32:10.000'
        DATE_AVG= '2010-07-23T15:03:07.000'
        DATE_END= '2010-07-23T15:34:04.000'
        TELESCOP= 'Hinode  '
        INSTRUME= 'EIS     '
        TARGET  = 'Active Region'
        STUD_ACR= 'HPW021_VEL_120x512v1'
        OBSTITLE= 'Scan of the core of AR 11089 as context for full CCD sit-and-stare.&'
        LINE_ID = 'Fe XII 195.119'
        MEASRMNT= 'intensity'
        BUNIT   = 'erg / (cm2 s sr)'
        CRVAL1  =   -430.9469060897827
        CRPIX1  =                    1
        CDELT1  =    1.996799945831299
        CTYPE1  = 'HPLN-TAN'
        CUNIT1  = 'arcsec  '
        CRVAL2  =   -693.4184875488281
        CRPIX2  =                    1
        CDELT2  =                  1.0
        CTYPE2  = 'HPLT-TAN'
        CUNIT2  = 'arcsec  '
        FOVX    =    119.8079967498779
        FOVY    =                512.0
        XCEN    =   -371.0429077148438
        YCEN    =   -437.4184875488281
        HGLN_OBS=                  0.0
        HGLT_OBS=    5.097799430847445
        DSUN_OBS=    151971313690.6726
        HISTORY fit using eispac 0.92.0 on 2021-09-27T13:58:29
        END
        """)
    return fits.Header.fromstring(raw_header, sep='\n')


@pytest.fixture
def test_eis_map(test_header):
    data = np.random.rand(test_header['naxis1'], test_header['naxis2'])
    return sunpy.map.Map(data, test_header)


def test_observatory(test_eis_map):
    assert test_eis_map.observatory == 'Hinode'


def test_wavelength(test_eis_map):
    assert test_eis_map.wavelength == 195.119*u.angstrom


def test_measurement(test_eis_map):
    assert test_eis_map.measurement == 'Fe XII 195.119 intensity'


def test_processing_level(test_eis_map):
    assert test_eis_map.processing_level == 3


def test_plot_settings(test_eis_map):
    assert test_eis_map.plot_settings['cmap'] == 'Blues_r'
    assert isinstance(test_eis_map.plot_settings['norm'], ImageNormalize)
    assert isinstance(test_eis_map.plot_settings['norm'].stretch, AsinhStretch)
    scale = test_eis_map.scale.axis2 / test_eis_map.scale.axis1
    assert test_eis_map.plot_settings['aspect'] == scale.decompose().value


def test_date(test_eis_map, test_header):
    # NOTE: these tests will change slightly with sunpy 3.1 decause of the existence of
    # date-beg, date-end, and date-avg keys
    # Case 1: date_obs and date_end exist
    assert test_eis_map.meta['date-beg'] == test_eis_map.meta['date_beg']
    assert test_eis_map.meta['date-end'] == test_eis_map.meta['date_end']
    assert test_eis_map.meta['date-obs'] == test_eis_map.meta['date_obs']
    assert test_eis_map.meta['date-avg'] == test_eis_map.meta['date_avg']
    assert test_eis_map.date.isot == test_eis_map.meta['date-avg']
    # Case 2: date_beg does not exist
    header = copy.deepcopy(test_header)
    del header['date_beg']
    new_map = sunpy.map.Map(test_eis_map.data, header)
    assert new_map.meta['date-beg'] == new_map.meta['date_obs']
    # Case 3: date_avg does not exist
    header = copy.deepcopy(test_header)
    del header['date_avg']
    new_map = sunpy.map.Map(test_eis_map.data, header)
    start = parse_time(new_map.meta['date-beg'], scale='utc')
    end = parse_time(new_map.meta['date-end'], scale='utc')
    avg = start + (end - start)/2
    assert new_map.meta['date-avg'] == avg.isot
    assert new_map.date == avg
    # Case 4: date_end and date_avg do not exist
    header = copy.deepcopy(test_header)
    del header['date_avg']
    del header['date_end']
    new_map = sunpy.map.Map(test_eis_map.data, header)
    assert new_map.date.isot == new_map.meta['date-obs']


def test_missing_date_raises_warning(test_eis_map, test_header):
    header = copy.deepcopy(test_header)
    del header['date_end']
    del header['date_obs']
    del header['date_beg']
    del header['date_avg']
    now = astropy.time.Time.now()
    new_map = sunpy.map.Map(test_eis_map.data, header)
    # This raises a slightly different warning in in sunpy>=3.1
    version = float('.'.join(sunpy.__version__.split('.')[:-1]))
    if version >= 3.1:
        from sunpy.util.exceptions import SunpyMetadataWarning
        expected_warning = SunpyMetadataWarning
    else:
        from sunpy.util.exceptions import SunpyUserWarning
        expected_warning = SunpyUserWarning
    with pytest.warns(expected_warning, match='Missing metadata for observation time'):
        assert new_map.date - now < 1*u.s

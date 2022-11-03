"""
Tests for `eispac.EISMap`
"""
import copy
from textwrap import dedent

import numpy as np
import pytest
from astropy.io import fits
import astropy.units as u
from astropy.visualization import ImageNormalize, AsinhStretch
import sunpy
import sunpy.map

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
    assert test_eis_map.measurement == 'intensity'


def test_nickname(test_eis_map):
    assert test_eis_map.nickname == 'Hinode EIS Fe XII 195.119'


def test_processing_level(test_eis_map):
    assert test_eis_map.processing_level == 3


def test_spatial_units(test_eis_map):
    assert test_eis_map.spatial_units[0] == u.arcsec
    assert test_eis_map.spatial_units[1] == u.arcsec


def test_waveunit(test_eis_map):
    assert test_eis_map.waveunit == u.Angstrom


@pytest.mark.parametrize('attr,key', [
    ('date_start', 'date_beg'),
    ('date_end', 'date_end'),
    ('date_average', 'date_avg'),
])
def test_date_props(test_eis_map, attr, key):
    assert getattr(test_eis_map, attr).isot == test_eis_map.meta[key]


def test_date(test_eis_map, test_header):
    # Case 1: date-average exists
    assert test_eis_map.date.isot == test_eis_map.meta['date_avg']
    # Case 2: date-average is None so default to date_obs
    header = copy.deepcopy(test_header)
    del header['date_beg']
    del header['date_end']
    del header['date_avg']
    new_map = sunpy.map.Map(test_eis_map.data, header)
    assert new_map.date.isot == new_map.meta['date_obs']
    # Case 3: date_avg is None so default to date_start
    header = copy.deepcopy(test_header)
    del header['date_avg']
    del header['date_end']
    del header['date_obs']
    new_map = sunpy.map.Map(test_eis_map.data, header)
    assert new_map.date.isot == new_map.meta['date_beg']
    # Case 4: date_end and date_avg do not exist
    header = copy.deepcopy(test_header)
    del header['date_avg']
    del header['date_beg']
    del header['date_obs']
    new_map = sunpy.map.Map(test_eis_map.data, header)
    assert new_map.date.isot == new_map.meta['date_end']


def test_plot_settings(test_eis_map):
    assert test_eis_map.plot_settings['cmap'] == 'Blues_r'
    assert isinstance(test_eis_map.plot_settings['norm'], ImageNormalize)
    assert isinstance(test_eis_map.plot_settings['norm'].stretch, AsinhStretch)
    scale = test_eis_map.scale.axis2 / test_eis_map.scale.axis1
    assert test_eis_map.plot_settings['aspect'] == scale.decompose().value

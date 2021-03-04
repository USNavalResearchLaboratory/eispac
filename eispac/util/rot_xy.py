
__all__ = ['rot_xy', 'test_rot_xy']

import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import Helioprojective, get_body_heliographic_stonyhurst
from sunpy.physics.differential_rotation import solar_rotate_coordinate
from sunpy.time import parse_time

def rot_xy(xcen, ycen, start_time, end_time):
    """
    Compute the rotation of a point on the sun

    xcen: solar x in arcsec (scalar)
    ycen: solar y in arcsec (scalar)
    start_time: starting time (will be parsed)
    end_time: end time (will be parsed)

    returns Tx, Ty in arcsec

    example

    > new = rot_xy(0, 0, start_time='01-JAN-2021 00:00', end_time='01-JAN-2021 01:00')
    > print(new.Tx, new.Ty)
    9.47188arcsec 0.0809565arcsec
    """

    start_time = parse_time(start_time)
    end_time = parse_time(end_time)

    c = SkyCoord(xcen*u.arcsec, ycen*u.arcsec, obstime=start_time,
                 observer='earth', frame=Helioprojective)

    observer = get_body_heliographic_stonyhurst('earth', end_time)

    new = solar_rotate_coordinate(c, observer=observer)

    return new

def test_rot_xy():
    
    start_time = '2012-09-24T10:50:26.000'
    xcen = -131.43
    ycen = -1.64
    end_time = '2012-09-24T11:50:26.000'    

    new = rot_xy(xcen, ycen, start_time, end_time)

    val = [-122.28026191424578, -1.7654771021022289]

    print(f' Computed {new.Tx.value:7.2f} {new.Ty.value:7.2f}')
    print(f' Expected {val[0]:7.2f} {val[1]:7.2f}')

    new = rot_xy(0, 0, start_time='2021-JAN-01 00:00', end_time='2021-JAN-01 01:00')
    print('\n ', new.Tx, new.Ty)    

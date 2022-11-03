"""
Convenience functions for rotating coordinates
"""
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import Helioprojective, get_body_heliographic_stonyhurst
from sunpy.physics.differential_rotation import solar_rotate_coordinate
from sunpy.time import parse_time

__all__ = ['rot_xy']


def rot_xy(xcen, ycen, start_time, end_time):
    """
    Compute the rotation of a point on the Sun using
    `~sunpy.physics.differential_rotation.solar_rotate_coordinate`

    Parameters
    ----------
    xcen: array-like
        Solar-x in arcsec (scalar)
    ycen: array-like
        Solar y in arcsec (scalar)
    start_time : any time format that can be parsed by `~sunpy.time.parse_time`
        Obstime of the coordinate represented by ``xcen`` and ``ycen``
    end_time: any time format that can be parsed by `~sunpy.time.parse_time`

    Returns
    -------
    new : `~astropy.coordinates.SkyCoord`
        Coordinate rotated by the time delta between ``start_time``
        and ``end_time``

    Examples
    --------
    >>> new = rot_xy(0, 0, start_time='2021-JAN-01 00:00', end_time='2021-JAN-01 01:00')
    >>> print(new.Tx, new.Ty)
    9.47188arcsec 0.0809565arcsec
    """
    start_time = parse_time(start_time)
    end_time = parse_time(end_time)
    c = SkyCoord(xcen*u.arcsec, ycen*u.arcsec, obstime=start_time,
                 observer='earth', frame=Helioprojective)
    observer = get_body_heliographic_stonyhurst('earth', end_time)
    new = solar_rotate_coordinate(c, observer=observer)
    return new    

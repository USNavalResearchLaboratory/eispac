from sunpy.net.dataretriever import GenericClient

from eispac.net.attrs import FileType

__all__ = ['EISClient']


class EISClient(GenericClient):
    """
    Provides access to the level 1 EIS data in HDF5 format.

    This data is hosted by the `Naval Research Laboratory <https://eis.nrl.navy.mil/>`__

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> import eispac.net
    >>> from eispac.net.attrs import FileType
    >>> results = Fido.search(a.Time('2020-11-09 00:00:00','2020-11-09 01:00:00'),
    ...                       a.Instrument('EIS'),
    ...                       a.Physobs.intensity,
    ...                       a.Source('Hinode'),
    ...                       a.Provider('NRL'),
    ...                       a.Level('1'))
    >>> results
    <sunpy.net.fido_factory.UnifiedResponse object at 0x7fe250ad3d00>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the EISClient:
           Start Time               End Time        Instrument  Physobs  Source Provider Level FileType
    ----------------------- ----------------------- ---------- --------- ------ -------- ----- --------
    2020-11-09 00:10:12.000 2020-11-09 00:10:12.999        EIS intensity HINODE      NRL     1     data
    2020-11-09 00:10:12.000 2020-11-09 00:10:12.999        EIS intensity HINODE      NRL     1     head
    <BLANKLINE>
    <BLANKLINE>
    >>> results = Fido.search(a.Time('2020-11-09 00:00:00','2020-11-09 01:00:00'),
    ...                       a.Instrument('EIS'),
    ...                       a.Physobs.intensity,
    ...                       a.Source('Hinode'),
    ...                       a.Provider('NRL'),
    ...                       a.Level('1'),
    ...                       FileType('head'))
    >>> results
    <sunpy.net.fido_factory.UnifiedResponse object at 0x7fd1baeeefd0>
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the EISClient:
           Start Time               End Time        Instrument  Physobs  Source Provider Level FileType
    ----------------------- ----------------------- ---------- --------- ------ -------- ----- --------
    2020-11-09 00:10:12.000 2020-11-09 00:10:12.999        EIS intensity HINODE      NRL     1     head
    <BLANKLINE>
    <BLANKLINE>
    """
    baseurl = r'https://eis.nrl.navy.mil/level1/hdf5/%Y/%m/%d/eis_%Y%m%d_%H%M%S.(\w){4}.h5'
    pattern = '{}/hdf5/{year:4d}/{month:2d}/{day:2d}/eis_{:8d}_{hour:2d}{minute:2d}{second:2d}.{FileType}.{}'
    
    @property
    def info_url(self):
        return 'https://eis.nrl.navy.mil/'
    
    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        return {
            attrs.Instrument: [('EIS', 'Extreme Ultraviolet Imaging Spectrometer')],
            attrs.Physobs: [('intensity', 'Spectrally resolved intensity in detector units')],
            attrs.Source: [('Hinode', 'The Hinode mission is a partnership between JAXA, NASA, and UKSA')],
            attrs.Provider: [('NRL', 'U.S. Naval Research Laboratory')],
            attrs.Level: [
                ('1', 'EIS: The EIS client can only return level 1 data. Level 0 EIS data is available from the VSO.')
            ],
            FileType: [('data', 'These files contain the actual intensity data'),
                       ('head', 'These files contain only the header metadata')],
        }
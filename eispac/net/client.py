from sunpy.net import attrs as a
from sunpy.net.dataretriever import GenericClient, QueryResponse
from sunpy.net.scraper import Scraper
from sunpy.time import TimeRange

from eispac.net.attrs import FileType

__all__ = ['EISClient']


class EISClient(GenericClient):
    """
    Provides access to the level 1 EIS data in HDF5 and FITS format.

    This data is hosted by the `Naval Research Laboratory <https://eis.nrl.navy.mil/>`__.

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
    ...                       a.Level('1'))  #doctest: +SKIP
    >>> results  #doctest: +SKIP
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the EISClient:
    Source: https://eis.nrl.navy.mil/
    <BLANKLINE>
           Start Time               End Time        ... Level   FileType
    ----------------------- ----------------------- ... ----- -----------
    2020-11-09 00:10:12.000 2020-11-09 00:10:12.999 ...     1   HDF5 data
    2020-11-09 00:10:12.000 2020-11-09 00:10:12.999 ...     1 HDF5 header
    2020-11-09 00:10:12.000 2020-11-09 00:10:12.999 ...     1        FITS
    <BLANKLINE>
    <BLANKLINE>
    >>> results = Fido.search(a.Time('2020-11-09 00:00:00','2020-11-09 01:00:00'),
    ...                       a.Instrument('EIS'),
    ...                       a.Physobs.intensity,
    ...                       a.Source('Hinode'),
    ...                       a.Provider('NRL'),
    ...                       a.Level('1'),
    ...                       FileType('HDF5 header'))  #doctest: +SKIP
    >>> results  #doctest: +SKIP
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the EISClient:
    Source: https://eis.nrl.navy.mil/
    <BLANKLINE>
           Start Time               End Time        ... Level   FileType
    ----------------------- ----------------------- ... ----- -----------
    2020-11-09 00:10:12.000 2020-11-09 00:10:12.999 ...     1 HDF5 header
    <BLANKLINE>
    <BLANKLINE>
    """
    baseurl_hdf5 = r'https://eis.nrl.navy.mil/level1/hdf5/%Y/%m/%d/eis_%Y%m%d_%H%M%S.(\w){4}.h5'
    pattern_hdf5 = '{}/{year:4d}/{month:2d}/{day:2d}/eis_{:8d}_{hour:2d}{minute:2d}{second:2d}.{FileType}'
    baseurl_fits = r'https://eis.nrl.navy.mil/level1/fits/%Y/%m/%d/eis_er_%Y%m%d_%H%M%S.fits'
    pattern_fits = '{}/{year:4d}/{month:2d}/{day:2d}/eis_er_{:8d}_{hour:2d}{minute:2d}{second:2d}.{FileType}'

    @property
    def info_url(self):
        return 'https://eis.nrl.navy.mil/'

    def search(self, *args, **kwargs):
        # NOTE: Search is overridden because URL and pattern depending on the filetype.
        # This enables multiple filetypes to be returned in the same query.
        metalist = []
        matchdict = self._get_match_dict(*args, **kwargs)
        all_filetypes = matchdict.get('FileType')
        for ft in all_filetypes:
            if 'h5' in ft:
                baseurl = self.baseurl_hdf5
                pattern = self.pattern_hdf5
            else:
                baseurl = self.baseurl_fits
                pattern = self.pattern_fits

            scraper = Scraper(baseurl, regex=True)
            tr = TimeRange(matchdict['Start Time'], matchdict['End Time'])
            filesmeta = scraper._extract_files_meta(tr, extractor=pattern, matcher={'FileType': ft})
            filesmeta = sorted(filesmeta, key=lambda k: k['url'])
            for i in filesmeta:
                rowdict = self.post_search_hook(i, matchdict)
                metalist.append(rowdict)

        return QueryResponse(metalist, client=self)

    def post_search_hook(self, i, matchdict):
        # This makes the final display names of the file types nicer
        filetype_mapping = {
            'data.h5': 'HDF5 data',
            'head.h5': 'HDF5 header',
            'fits': 'FITS',
        }
        rd = super().post_search_hook(i, matchdict)
        rd['FileType'] = filetype_mapping[rd['FileType']]
        return rd

    @classmethod
    def register_values(cls):
        return {
            a.Instrument: [('EIS', 'Extreme Ultraviolet Imaging Spectrometer')],
            a.Physobs: [('intensity', 'Spectrally resolved intensity in detector units')],
            a.Source: [('Hinode', 'The Hinode mission is a partnership between JAXA, NASA, and UKSA')],
            a.Provider: [('NRL', 'U.S. Naval Research Laboratory')],
            a.Level: [
                ('1', 'EIS: The EIS client can only return level 1 data. Level 0 EIS data is available from the VSO.')
            ],
            FileType: [('data.h5', 'These files contain the actual intensity data in HDF5 format.'),
                       ('head.h5', 'These files contain only the header metadata in HDF5 format.'),
                       ('fits', 'These files contain both data and metadata in FITS format')],
        }

    @classmethod
    def _attrs_module(cls):
        # Register EIS specific attributes with Fido
        return 'eispac', 'eispac.net.attrs'

    @classmethod
    def _can_handle_query(cls, *query):
        """
        Check if this client can handle a given Fido query.
        Returns
        -------
        bool
            True if this client can handle the given query.
        """
        required = {a.Time, a.Instrument, a.Source}
        optional = {a.Provider, a.Physobs, a.Level, FileType}
        return cls.check_attr_types_in_query(query, required, optional)

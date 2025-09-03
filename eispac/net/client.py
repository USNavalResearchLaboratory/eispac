__all__ = ['EISClient']

import copy
from pathlib import Path
import numpy as np
from sunpy import __version__ as sunpy_ver
from sunpy.net import attrs as a
from sunpy.net.attr import SimpleAttr
from sunpy.net.dataretriever import GenericClient, QueryResponse
from sunpy.net.base_client import QueryResponseRow
from sunpy.net.scraper import Scraper
from sunpy.time import TimeRange
from sunpy.util.parfive_helpers import Downloader

from eispac.net.attrs import FileType

class EISClient(GenericClient):
    """Provides access to the level 1 EIS data in HDF5 and FITS format.

    The EIS level-1 HDF5 files come in pairs of "data.h5" and "head.h5" files
    containing, respectively, the EUV spectra and header/metadata for a single
    observation. By default, this client only returns entries for the data files
    and fetching a data.h5 file will automatically download the matching header.
    The original level-1 FITS files are also available, should a user wish to 
    analyze the data using older IDL routines. Fetching a FITS entry will download 
    both the "l1" (data) and "er" (errors) files produced by the "EIS_PREP" 
    routine in SolarSoft/IDL.

    This data is processed and hosted by the U.S. Naval Research Laboratory (NRL)
    A mirror of the HDF5 files is hosted by the Mullard Space Science Laboratory
    (MSSL) in the UK.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> import eispac.net
    >>> results = Fido.search(a.Time('2020-11-09 00:00:00','2020-11-09 01:00:00'),
    ...                       a.Instrument('EIS'),
    ...                       a.Source('Hinode'),
    ...                       a.Provider('NRL'))  #doctest: +SKIP
    >>> results  #doctest: +SKIP
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the EISClient:
    Source: https://eis.nrl.navy.mil/
    <BLANKLINE>
           Start Time               End Time        ... Level   FileType
    ----------------------- ----------------------- ... ----- -----------
    2020-11-09 00:10:12.000 2020-11-09 00:10:12.999 ...     1   HDF5 data
    <BLANKLINE>
    <BLANKLINE>
    >>> results = Fido.search(a.Time('2020-11-09 00:00:00','2020-11-09 01:00:00'),
    ...                       a.Instrument('EIS'),
    ...                       a.Source('Hinode'),
    ...                       a.Provider('NRL'),
    ...                       a.eispac.FileType('Any'))  #doctest: +SKIP
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
    """

    current_provider = 'nrl'
    nrl_topurl = r'https://eis.nrl.navy.mil/level1/'
    mssl_topurl = r'https://vsolar.mssl.ucl.ac.uk/eispac/'
    baseurl_hdf5 = r'hdf5/%Y/%m/%d/eis_%Y%m%d_%H%M%S.(\w){4}.h5'
    pattern_hdf5 = '{}/{year:4d}/{month:2d}/{day:2d}/eis_{:8d}_{hour:2d}{minute:2d}{second:2d}.{FileType}'
    format_hdf5 = r'hdf5/{{year:4d}}/{{month:2d}}/{{day:2d}}/'+ \
                  r'eis_{{:8d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}.{{FileType}}'
    baseurl_fits = r'fits/%Y/%m/%d/eis_l1_%Y%m%d_%H%M%S.fits'
    pattern_fits = '{}/{year:4d}/{month:2d}/{day:2d}/eis_l1_{:8d}_{hour:2d}{minute:2d}{second:2d}.{FileType}'
    format_fits = r'fits/{{year:4d}}/{{month:2d}}/{{day:2d}}/'+ \
                  r'eis_l1_{{:8d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}.{{FileType}}'

    @property
    def info_url(self):
        if self.current_provider.lower() == 'nrl':
            return 'https://eis.nrl.navy.mil/'
        if self.current_provider.lower() == 'mssl':
            return 'https://vsolar.mssl.ucl.ac.uk/eispac/hdf5/'

    def search(self, *args, **kwargs):
        # NOTE: Search is overridden because URL and pattern depend on filetype.
        # This enables multiple filetypes to be returned in the same query.
        metalist = []
        matchdict = self._get_match_dict(*args, **kwargs)
        t_range = TimeRange(matchdict['Start Time'], matchdict['End Time'])
        all_filetypes = matchdict.get('FileType')
        if len(all_filetypes) > 1:
            # When FileType is NOT specified, default to 'data.h5' (first index)
            all_filetypes = [all_filetypes[0]]
        if all_filetypes[0].lower() == 'any':
            all_filetypes = ['data.h5', 'head.h5', 'fits']
        matchdict['FileType'] = all_filetypes # copy changes back to matchdict
        for f_type in all_filetypes:
            if matchdict.get('Provider')[0].lower() == 'nrl':
                self.current_provider = 'nrl'
                topurl = self.nrl_topurl
            elif matchdict.get('Provider')[0].lower() == 'mssl':
                self.current_provider = 'mssl'
                if 'fits' in f_type:
                    print('NOTICE: level-1 FITS files are only available from NRL')
                    continue
                else:
                    topurl = self.mssl_topurl

            if 'h5' in f_type:
                baseurl = topurl + self.baseurl_hdf5
                pattern = self.pattern_hdf5
                parse_format = topurl + self.format_hdf5
            elif 'fits' in f_type:
                baseurl = topurl + self.baseurl_fits
                pattern = self.pattern_fits
                parse_format = topurl + self.format_fits
            else:
                continue # skip unknown filetypes

            if sunpy_ver < '6.1':
                # Old regex Scraper
                scraper = Scraper(baseurl, regex=True)
                filesmeta = scraper._extract_files_meta(t_range, extractor=pattern, 
                                                        matcher={'FileType': f_type})
            else:
                # New parse Scraper (added in SunPy 6.1 on 2025-02-24)
                scraper = Scraper(format=parse_format)
                filesmeta = scraper._extract_files_meta(t_range, 
                                                        matcher={'FileType': f_type})
            filesmeta = sorted(filesmeta, key=lambda k: k['url'])
            for i in filesmeta:
                rowdict = self.post_search_hook(i, matchdict)
                metalist.append(rowdict)

        if len(metalist) > 0:
            return QueryResponse(metalist, client=self)
        else:
            return QueryResponse(None, client=self)
    
    def fetch(self, qres, path=None, overwrite=False,
              progress=True, downloader=None, wait=True, **kwargs):
        """
        Download a set of results.

        Parameters
        ----------
        qres : `~sunpy.net.dataretriever.QueryResponse`
            Results to download.
        path : `str` or `pathlib.Path`, optional
            Path to the download directory, or file template including the
            ``{file}`` string which will be replaced with the filename.
        overwrite : `bool` or `str`, optional
            Determine how to handle downloading if a file already exists with the
            same name. If `False` the file download will be skipped and the path
            returned to the existing file, if `True` the file will be downloaded
            and the existing file will be overwritten, if ``'unique'`` the filename
            will be modified to be unique.
        progress : `bool`, optional
            If `True` show a progress bar showing how many of the total files
            have been downloaded. If `False`, no progress bar will be shown.
        downloader : `parfive.Downloader`, optional
            The download manager to use.
        wait : `bool`, optional
            If `False` ``downloader.download()`` will not be called. Only has
            any effect if ``downloader`` is not `None`.

        Returns
        -------
        `parfive.Results`
        """
        # NOTE: fetch is overwritten since we often need to download pairs of
        #       files for a single obs (data.h5 & head.h5 or l1 & er FITS) 

        if path is not None:
            path = Path(path)

        if isinstance(qres, QueryResponseRow):
            qres = qres.as_table()

        urls = []
        if len(qres):
            urls = list(qres['url'])

        filenames = [url.split('/')[-1] for url in urls]

        paths = self._get_full_filenames(qres, filenames, path)

        # If downloading data.h5 files, always download the head.h5 files too
        data_h5_in_list = any(['data.h5' in NAME for NAME in filenames])
        head_h5_in_list = any(['head.h5' in NAME for NAME in filenames])
        if data_h5_in_list and not head_h5_in_list:
            # Add missing head.h5 files to the list
            urls, paths = self._add_missing_file_pairs(urls, paths, 'data.h5', 'head.h5')

        # If downloading FITS files, always download both the l1 and er files
        fits_in_list = any(['eis_l1_' in NAME for NAME in filenames])
        if fits_in_list:
            # Add missing "eis_er_*.fits" files to the list
            urls, paths = self._add_missing_file_pairs(urls, paths, 'eis_l1_', 'eis_er_')
        
        dl_set = True
        if not downloader:
            dl_set = False
            downloader = Downloader(progress=progress, overwrite=overwrite)

        for url, filename in zip(urls, paths):
            downloader.enqueue_file(url, filename=filename, **self.enqueue_file_kwargs)

        if dl_set and not wait:
            return

        return downloader.download()

    def _add_missing_file_pairs(self, urls, paths, find_str, replace_str):
        # Appending to and sorting file lists to include missing pairs 
        new_urls = copy.deepcopy(urls)
        new_paths = copy.deepcopy(paths)
        for i in range(len(urls)):
            new_urls.append(np.str_(urls[i].replace(find_str, replace_str)))
            new_paths.append(Path(str(paths[i]).replace(find_str, replace_str)))

        urls = sorted(new_urls)
        paths = sorted(new_paths)

        return urls, paths

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
            a.Provider: [('NRL', 'U.S. Naval Research Laboratory'),
                         ('MSSL', 'Mullard Space Science Laboratory (UK)')],
            a.Level: [
                ('1', 'EIS: processed data with corrected metadata. Use the VSO for Level-0 data.')
            ],
            FileType: [('data.h5', 'Level-1 EUV spectra in HDF5 format'),
                       ('head.h5', 'Level-1 header & metadata in HDF5 format'),
                       ('fits', 'Pairs of "l1" and "er" files in FITS format'),
                       ('any', 'Retrieve any/all file types available')],
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
        
        # Check that all query attrs are accepted by this client
        if not cls.check_attr_types_in_query(query, required, optional):
            return False

        # Validate the actual attribute values
        regattrs_dict = cls.register_values()
        for key in regattrs_dict:
            all_vals = [i[0].lower() for i in regattrs_dict[key]]
            for x in query:
                if (isinstance(x, key) and issubclass(key, SimpleAttr) 
                and str(x.value).lower() not in all_vals):
                    return False
        return True

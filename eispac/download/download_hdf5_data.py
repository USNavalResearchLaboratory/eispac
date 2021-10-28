#!/usr/bin/env python

import os
import parfive
# import urllib.request, requests, wget


class download_hdf5_data:
    """
    An object for downloading EIS HDF5 level1 files

    input is parsed to construct the remote and local filenames, curl is spawned to download
    the files (what if it doesn't exist?), if files exist locally they are skipped.

    Parameters
    ----------
    filename :  `str` or `list`
        An EIS filename. Such as, eis_l0_20200311_213413.fits, eis_l1_20200311_213413.fits.gz,
        /some_path/eis_l0_20200311_213413.fits. Can be a single path or a list of paths
    datetree : `bool`
        Create a local path organized by date (YYYY/MM/DD)
    local_top : `str`
        Top of the local path (e.g., data_eis)
    nodata : `bool`
        Don't download data files
    nohead : `bool`
        Don't download head files
    overwrite : `bool`
        Download even if file exists locally
    headonly : `bool`
        Equivalent to ``nodata`` + ``overwrite``
    max_conn : int
        Max number of download connections that parfive will use

    """

    def __init__(self, filename=None, local_top='data_eis', datetree=False,
                 nodata=False, nohead=False, overwrite=False, headonly=False,
                 max_conn=2):
        self.top_url = 'https://eis.nrl.navy.mil/level1/hdf5/'
        self.local_top = local_top
        self.datetree = datetree
        self.nodata = nodata
        self.nohead = nohead
        self.overwrite = overwrite
        self.max_conn = max_conn
        if headonly:
            self.nodata = True
            self.overwrite = True
        if local_top is None:
            self.local_top = os.getcwd()
        elif local_top.casefold() == 'cwd':
            self.local_top = os.getcwd()

        if filename is not None:
            self.process_input(filename)

    def process_input(self, this_input):
        """
        input can be a filename or a list of filenames
        """
        # make input into a list so we can iterate over it
        if not isinstance(this_input, list):
            this_input = [this_input]
        # process input and download
        for this_file in this_input:
            self.construct_filenames(this_file)
            self.download()

    def construct_filenames(self, input_filename):
        """
        convert input into remote and local filenames
        """
        f, year, month, day = self.eis_filename = self.parse_input_filename(input_filename)
        self.f_data = f + '.data.h5'
        self.f_head = f + '.head.h5'
        self.date_path = os.path.join(year, month, day)
        self.remote_data = os.path.join(self.top_url, self.date_path, self.f_data)
        self.remote_head = os.path.join(self.top_url, self.date_path, self.f_head)
        self.remote_data = self.remote_data.replace('\\', '/') # ensure unix format
        self.remote_head = self.remote_head.replace('\\', '/')
        self.local_data = self.f_data
        self.local_head = self.f_head
        if self.datetree:
            self.local_path = os.path.join(self.local_top, self.date_path)
        else:
            self.local_path = self.local_top
        self.local_data = os.path.join(self.local_path, self.local_data)
        self.local_head = os.path.join(self.local_path, self.local_head)

    def parse_input_filename(self, input_filename):
        """
        convert an eis fits filename into an hdf5 filename; extract year, month, day
        """
        # Extract and format the filename
        # For exmaple, /path/eis_l0_20200311_213413.fits.gz
        f = os.path.basename(input_filename)
        f = f.split('.')[0] # eis_l0_20200311_213413
        f = f.replace('_l0_', '_') # eis_20200311_213413
        f = f.replace('_l1_', '_') # eis_20200311_213413
        # Extract year, month, day from filename
        d = f.split('_')[1] # 20200311
        year = d[0:4]
        month = d[4:6]
        day = d[6:8]
        return f, year, month, day

    def download(self):
        """
        download data and head files, unless told not to
        """
        if not self.nodata:
            self.download_file(self.remote_data, self.local_path, self.f_data)
        if not self.nohead:
            self.download_file(self.remote_head, self.local_path, self.f_head)

    def download_file(self, remote_filepath, local_dir, local_name):
        """
        Use parfive to download the file
        """
        local_filepath = os.path.join(local_dir, local_name)
        # check that local directory exists and create as needed
        self.check_local_dir(local_filepath)
        if self.overwrite or not os.path.isfile(local_filepath):
            try:
                print(f'+ downloading {remote_filepath} -> {local_filepath}')
                dl = parfive.Downloader(max_conn=self.max_conn, overwrite=True)
                dl.enqueue_file(remote_filepath, path=local_dir,
                                filename=local_name+'.part', verify_ssl=False)
                dl_result = dl.download()
                if len(dl_result.errors) == 0:
                    os.replace(local_filepath+'.part', local_filepath)
                else:
                    print(f' ERROR: Incomplete file download.')
                print('')
            except:
                print(f' ! Error trying to download: {remote_filepath}')
        else:
            print(f' + {local_filepath} exists, skipping download')

    def check_local_dir(self, local_file):
        """
        check if local dir exits. If not, create it
        """
        local_dir = os.path.dirname(local_file)
        if local_dir != '':
            if not os.path.isdir(local_dir):
                os.makedirs(local_dir)
                print(' + created ' + local_dir)

if __name__ == '__main__':

    # do nothing
    o = download_hdf5_data(filename='eis_l0_20200311_213413.fits')

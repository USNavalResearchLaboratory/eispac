#!/usr/bin/env python

import glob
import os
import numpy as np
import h5py
import pathlib
import eispac
from eispac.core.read_wininfo import read_wininfo as eis_read_wininfo

class EISTemplateLocator:
    """Helper class used by the "eis_browse_templates" GUI

    Finds fit templates matching all data windows in a given EIS observation
    and generates lists of the filepaths. Not intended for direct use.
    For a more user-friendly function with similar capabilities, see
    `~eispac.core.match_templates`

    Parameters
    ----------
    filename_head : str, optional
        Filepath to a EIS Level-1 HDF5 file. Both data and header filepaths are
        accepted.

    verbose :bool, optional
        If set to True, will automatically print a list of all found templates.
        Default is False.

    ignore_local : bool, optional
        If set to True, will look for templates in either SSW or those
        distributed with EISPAC. If set to False, will instead ONLY look for
        templates in a local directory named either "eis_template_dir" or
        "eis_fit_templates". Default is True.
    """

    def __init__(self, filename_head=None, verbose=False, ignore_local=True):
        self.filename_head = None
        self.template_list = None
        self.n_template_files = 0
        self.wininfo = None
        self.EIS_TEMPLATE_DIR = None
        self.text_list = None
        self.ignore_local = ignore_local
        self.delta_wave = 0.25

        if filename_head is not None:
            self.filename_head = self.parse_input_filename(filename_head)

            self.find_templates()
            self.read_wininfo()
            self.match_templates()
            self.make_text_list()
            if verbose:
                self.print_text_list()

    def parse_input_filename(self, input_filename):
        """Read an EIS level-1 filename and determine the header filepath"""
        path = os.path.dirname(input_filename)
        base = os.path.basename(input_filename).replace('.data.h5','.head.h5')
        return os.path.join(path, base)

    def find_templates(self):
        """Find all fit template files available"""
        if os.path.isdir('eis_template_dir') and (not self.ignore_local):
            # look for local templates (old default, before 2021-10-14)
            self.EIS_TEMPLATE_DIR = 'eis_template_dir'
        elif os.path.isdir('eis_fit_templates') and (not self.ignore_local):
            # look for local templates (current default, after 2021-10-14)
            self.EIS_TEMPLATE_DIR = 'eis_fit_templates'
        elif os.getenv('EIS_TEMPLATE_DIR') is not None:
            # look for an environment variable
            self.EIS_TEMPLATE_DIR = os.getenv('EIS_TEMPLATE_DIR')
        else:
            # look for the templates distributed with the package
            # root = str(pathlib.Path(__file__).parent.parent.absolute())
            root = os.path.dirname(eispac.__file__)
            path = os.path.join(root, 'data', 'templates')
            self.EIS_TEMPLATE_DIR = path

        # find all of the template files
        search = os.path.join(self.EIS_TEMPLATE_DIR, '*.template.h5')
        self.template_files = glob.glob(search)
        self.n_template_files = len(self.template_files)

        if self.n_template_files > 0:
            self.template_files = sorted(self.template_files)
        else:
            print(' ! no template files found')

    def read_wininfo(self):
        """"Read the window information from the header file"""
        if os.path.isfile(self.filename_head):
            self.wininfo = eis_read_wininfo(self.filename_head)
            for iwin in range(len(self.wininfo)):
                self.wininfo[iwin].wvl_min -= self.delta_wave
                self.wininfo[iwin].wvl_max += self.delta_wave
            if self.wininfo is None:
                print(' ! wininfo not read')
        else:
            print(' ! EIS head file not read')

    def match_templates(self):
        """Match each observation window and with all relevant fit templates"""
        if (self.n_template_files == 0) or (self.wininfo is None):
            return

        template_list = []
        for tf in self.template_files:
            with h5py.File(tf, 'r') as f:
                wmin = f['/template/wmin'][0]
                wmax = f['/template/wmax'][0]
                match, = np.where((wmin >= self.wininfo.wvl_min) & (wmax <= self.wininfo.wvl_max))
                if len(match) > 0:
                    for m in match:
                        template_list.append((tf, m))

        # template list is (template file, iwin)
        template_list = sorted(template_list, key=lambda x: x[1])
        self.template_list = template_list

    def make_text_list(self):
        """Generate a text list of all templates found"""
        if (self.n_template_files == 0) or (self.wininfo is None):
            return

        txt = []
        wininfo = self.wininfo
        for t in self.template_list:
            f = os.path.basename(t[0])
            iwin = t[1]
            txt.append(f'{f}{iwin:6d}{wininfo[t[1]].wvl_min:8.2f}{wininfo[t[1]].wvl_max:8.2f}')
        self.text_list = txt

    def print_text_list(self):
        """Print the list of found templates"""
        for t in self.text_list:
            print(t)

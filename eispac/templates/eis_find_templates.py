#!/usr/bin/env python

import glob
import os
import numpy as np
import h5py
from pathlib import Path
import eispac
from eispac.core.read_wininfo import read_wininfo as eis_read_wininfo

class eis_find_templates:
    def __init__(self, filename_head=None, verbose=False, ignore_local=False):
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
            if verbose: self.print_text_list()

    def parse_input_filename(self, input_filename):
        path = os.path.dirname(input_filename)
        base = os.path.basename(input_filename).replace('.data.h5','.head.h5')
        return os.path.join(path, base)

    def find_templates(self):
        # find the templates
        if os.path.isdir('eis_template_dir') and (not self.ignore_local):
            # look for local templates
            self.EIS_TEMPLATE_DIR = 'eis_template_dir'
        elif os.getenv('EIS_TEMPLATE_DIR') is not None:
            # look for an environment variable
            self.EIS_TEMPLATE_DIR = os.getenv('EIS_TEMPLATE_DIR')
        else:
            # look for the templates distributed with the package
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
        # read the window information from the file
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
        # match windows and templates
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
        for t in self.text_list:
            print(t)

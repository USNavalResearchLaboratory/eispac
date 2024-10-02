#!/usr/bin/env python
"""EIS as-run catalog search PyQt5 gui.

A python gui to search the EIS as-run catalog. The intention is to
be somewhat like the IDL routine eis_cat.pro, but with fewer
features, at least initially. Very much a work in progress.

(2017-Apr-12) First working version added to git.
(2020-Dec-16) Now tries to find SSW if the 'SSW' environment variable is missing
(2020-Dec-18) If no database is found, ask the user if they want to download it.
(2021-Oct-22) Major update adding in more search criteria and faster queries.
(2023-Jul-11) Added data mirrors for eis_cat.sqlite and HDF5 files
(2023-Sep-23) Added viewing context images from MSSL
(2024-May-07) Removed PyQt4 support, cleaned-up code and fixed a search bug
"""
__all__ = ['eis_catalog']

import sys
import time
import os
import re
import urllib
import ssl
import certifi
import sqlite3
from datetime import datetime, timedelta
import numpy as np
from PyQt5 import QtCore, QtWidgets, QtGui, QtNetwork
from eispac.download import eis_obs_struct
from eispac.download.download_hdf5_data import download_hdf5_data
from eispac.download.download_db import download_db

def get_remote_image_dir(filename):
    """Parse a Level-0 filename and get the remote image dir"""
    # Note: level-0 files have the form, eis_l0_YYYYMMDD_hhmmss.fits
    date_str = filename.split('_')[2]
    year_str = date_str[0:4]
    month_str = date_str[4:6]
    day_str = date_str[6:8]
    file_dir = year_str+'/'+month_str+'/'+day_str+'/'+filename+'/'

    # Assemble and return the URL
    base_url = 'https://solarb.mssl.ucl.ac.uk/SolarB/DEV/eis_gifs/'
    return base_url + file_dir

class Top(QtWidgets.QWidget):

    def __init__(self, dbfile, parent=None):
        super(Top, self).__init__(parent)
        self.file_list = None
        self.selected_file = None
        self.selected_info = []
        self.default_filename = 'eis_filelist.txt'
        self.default_start_time = '2018-05-29 00:00' # '29-May-2018 00:00'
        self.default_end_time = '2018-05-29 23:59' # '29-May-2018 23:59'
        self.default_button_width = 150 #165 #130
        self.default_topdir = os.path.join(os.getcwd(), 'data_eis')
        self.dbfile = dbfile
        self.db_loaded = False
        self.context_imgNX = 768 #512
        self.context_imgNY = 768 #512

        # Font settings
        self.default_font = QtGui.QFont()
        self.small_font = QtGui.QFont()
        self.info_detail_font = QtGui.QFont("Courier New", 9)
        self.default_font.setPointSize(11)
        self.small_font.setPointSize(9)

        # Dict of search criteria to include as options in the drop-down list
        # Note: the keys:value pairs give the mapping of gui_label:sqlite_col
        self.criteria = {'Date Only':'date_obs', 'Study ID':'study_id',
                         'Study Acronym':'stud_acr', 'HOP ID':'jop_id',
                         'Obs. Title':'obstitle',
                         'Raster ID':'rast_id', 'Raster Acr.':'rast_acr',
                         'Triggered Obs':None,
                         'Target':'target', 'Science Obj.':'sci_obj',
                         'Timeline ID':'tl_id'}
        # Filter drop-down lists. given as gui_label:filter_value pairs
        self.rast_types = {'Any':None, 'Scan (0)':0, 'Sit-and-Stare (1)':1}
        self.slit_slot = {'Any':None, 'Slit only (0 & 2)':[0,2],
                          '1" slit (0)':[0], '2" slit (2)':[2],
                          'Slot only (3 & 1)':[3,1],
                          '40" slot (3)':[3], '266" slot (1)':[1]}

        self.rtype_dict = {0:'0 (scan)', 1:'1 (sit-and-stare)'}
        self.sindex_dict = {0:'0 (1" slit)', 2:'2 (2" slit)',
                            3:'3 (40" slot)', 1:'1 (266" slot)'}
        
        # Check for EIS database
        if os.path.isfile(self.dbfile):
            self.d = eis_obs_struct.EIS_DB(self.dbfile)
            self.db_loaded = True
        else:
            # Ask if the user wants to download a copy of the database
            ask_db = QtWidgets.QMessageBox.question(self, 'Missing database',
                        'No EIS as-run database found\n'
                        'Would you like to download a local copy?',
                        QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
            if ask_db == QtWidgets.QMessageBox.Yes:
                self.dbfile = download_db()
                if not str(self.dbfile).endswith('.part'):
                    self.d = eis_obs_struct.EIS_DB(self.dbfile)
                else:
                    print('Failed to download EIS database!')
            else:
                print('ERROR: No EIS as-run database found!')

        # Initialze the network manager for downloading images
        self.manager = QtNetwork.QNetworkAccessManager()
        self.manager.finished.connect(self.on_finished)

        self.init_ui()

    def init_ui(self):
        """Manage everything."""
        self.grid = QtWidgets.QGridLayout(self)
        self.gui_row = 0
        # self.setStyleSheet("QLabel{font-size: 11pt;}")

        # Quit & Update DB buttons (with filename & timestamp)
        self.top_menu() # very top

        # Selecting search criteria
        self.select_dates() # top left
        self.select_primary() # top center

        # Filter values
        self.set_filters() # top center

        # Catalog info
        self.catalog_table() # middle left

        # Info for a single search result
        self.details() # right

        # Bottom stuff
        self.save_options() # bottom left

        # And away we go
        self.setLayout(self.grid)
        self.setGeometry(50, 100, 1800, 800) #(50, 100, 1800, 800)
        self.setWindowTitle('EIS As-Run Catalog Information')
        self.event_help()
        self.show()

    def event_quit(self):
        QtWidgets.QApplication.instance().quit()

    def top_menu(self):
        """Basic menu option."""
        self.quit = QtWidgets.QPushButton('Quit')
        self.quit.setFixedWidth(self.default_button_width)
        self.quit.setFont(self.default_font)
        self.quit.clicked.connect(self.event_quit)

        self.download_db = QtWidgets.QPushButton('Update Database', self)
        self.download_db.setFixedWidth(self.default_button_width)
        self.download_db.setFont(self.default_font)
        self.download_db.clicked.connect(self.event_download_db)

        self.db_source_box = QtWidgets.QComboBox()
        self.db_source_box.addItems(['NASA (source)', 'NRL (mirror)'])
        self.db_source_box.setFixedWidth(self.default_button_width) # or 150
        self.db_source_box.setFont(self.default_font)

        self.db_info = QtWidgets.QLabel(self)
        self.db_info.setFixedWidth(4*self.default_button_width)
        self.db_info.setFont(self.small_font)
        if os.path.isfile(self.dbfile):
            self.update_db_file_label()
        else:
            self.db_info.setText('Unable to locate DB: ' + self.dbfile)

        self.grid.addWidget(self.quit, self.gui_row, 0)
        self.grid.addWidget(self.download_db, self.gui_row, 1)
        self.grid.addWidget(self.db_source_box, self.gui_row, 2)
        self.grid.addWidget(self.db_info, self.gui_row, 3, 1, 4)

        self.gui_row += 1

    def update_db_file_label(self):
         t = 'DB file: '+os.path.abspath(self.dbfile)+'\nDownload date: '+ \
                time.ctime(os.path.getmtime(self.dbfile)).lstrip().rstrip()
         self.db_info.setText(t)

    def event_download_db(self):
        self.tabs.setCurrentIndex(0) #switch to details tab
        self.info_detail.clear()
        self.table_m.clearContents()
        self.table_m.setRowCount(1)
        db_source = self.db_source_box.currentText()
        if db_source.lower().startswith('nrl'):
            db_remote_text = f'https://eis.nrl.navy.mil/level1/db/eis_cat.sqlite'
        else:
            db_remote_text = f'https://hesperia.gsfc.nasa.gov/ssw/hinode/eis' \
                            +f'/database/catalog/eis_cat.sqlite'
        info = (f'Downloading eis_cat.sqlite.\n'
                f'   Remote: {db_remote_text}\n'
                f'   Local: {os.path.abspath(self.dbfile)}\n\n'
                f'Please wait...')
        self.info_detail.append(info)
        QtWidgets.QApplication.processEvents() # update gui while user waits
        if self.db_loaded:
            self.d.cur.close()
            self.d.conn.close()
        self.d = 0
        self.dbfile = download_db(os.path.dirname(self.dbfile), source=db_source)
        if self.dbfile.endswith('.part'):
            self.info_detail.append('\nERROR: Database download failed! '
                                    +'Please check your internet connection '
                                    +' or try a different source.')
        else:
            self.d = eis_obs_struct.EIS_DB(self.dbfile)
            self.update_db_file_label()
            self.info_detail.append('\nComplete')
            self.db_loaded = True

    def event_help(self):
        """Put help info in details window."""
        self.info_help.clear()
        help_text = """EIS As-Run Catalog Search Tool

### Searching for Obervations

* Please use ISO format dates (YYYY-MM-DD HH:MM). If only the
  start date is provided, the end will be set for 24 hours later.
  WARNING: "Date Only" searches over the entire mission can take
  a VERY long time, please be patient.

* Primary search criteria descriptions (with catalog column names):

  Study ID (study_id)
      ID number for specific observation plan (line list, raster
      steps, etc.). Studies may repeated throughout the mission.

  Study Acronym (stud_acr)
      Short text label for the study. You only need to input a
      few characters (case is ignored).

  HOP ID (jop_id)
      Hinode Operation Plan ID. Assigned to observations that were
      coordinated with other telescopes or spacecraft missions.
      Also known as "JOP ID" (Joint Obs. Program ID)

  Obs. Title (obstitle)
      Observation title. Just a word or two is enough, results will
      be shown for all titles containing the input text.

  Raster ID (rast_id) and Raster Acr. (rast_acr)
      ID number and short text label for the raster program used in
      the study. Rasters may be shared by multiple studies.
      
  Triggered Obs
      Any raster run in response to an on-board trigger. There 
      are three triggers that can be run: XRT flare (tl_id=1),
      EIS flare (tl_id=3), & EIS bright point (tl_id=4). 
      Selecting this the same as using Timeline ID = 1, 3, 4 

  Target (target)
      Main observation target (e.g., Active Region, Quiet Sun).

  Science Obj. (sci_obj)
      Target phenomena (e.g., AR, QS, BP, LMB)
  
  Timeline ID (tl_id)
      ID number for a unique set of contiguous observations. This
      could be a single raster or a set of multiple rasters run in
      sequence. Timeline IDs of 1, 3, & 4 are special (see above), 
      all other values should be unique over the life of EIS.

  Please Note: If "Target" and "Science Obj." are not defined by
  the EIS planner, default values "Quiet Sun" & "QS" are assigned,
  regardless of the actual observation target.

  * Multiple ID numbers may be searched at the same time by separating
    the numbers with "," (e.g. "1, 3, 4"). Ranges of numbers can be 
    searched by separating the start & end values with "-" (e.g. "1-4")

### Downloading files

* If the "Use Date Tree" box is checked, files will be downloaded
  into subdirectories organized by date (../YYYY/MM/DD/)
"""
        self.info_help.append(help_text)
        self.info_help.verticalScrollBar().setValue(
                                self.info_help.verticalScrollBar().minimum())

    def select_dates(self):
        """Set time range and make button for running the search"""
        title = QtWidgets.QLabel(self)
        title.setText('Select Time Range')
        title.setFont(self.default_font)
        self.grid.addWidget(title, self.gui_row, 0, 1, 2)
        self.gui_row += 1

        start_t = QtWidgets.QLabel(self)
        start_t.setText('Start Time')
        start_t.setFont(self.default_font)
        self.start_time = QtWidgets.QLineEdit(self)
        self.start_time.setFixedWidth(self.default_button_width)
        self.start_time.setText(self.default_start_time)
        self.start_time.setFont(self.default_font)
        self.grid.addWidget(start_t, self.gui_row, 0)
        self.grid.addWidget(self.start_time, self.gui_row, 1)
        self.gui_row += 1

        end_t = QtWidgets.QLabel(self)
        end_t.setText('End Time')
        end_t.setFont(self.default_font)
        self.end_time = QtWidgets.QLineEdit(self)
        self.end_time.setFixedWidth(self.default_button_width)
        self.end_time.setText(self.default_end_time)
        self.end_time.setFont(self.default_font)
        self.grid.addWidget(end_t, self.gui_row, 0)
        self.grid.addWidget(self.end_time, self.gui_row, 1)
        self.gui_row += 1

        time_recent = QtWidgets.QPushButton('Last 3 Weeks', self)
        time_recent.setFixedWidth(self.default_button_width)
        time_recent.setFont(self.default_font)
        self.grid.addWidget(time_recent, self.gui_row, 0, 1, 6)
        time_recent.clicked.connect(self.event_time_recent)

        time_mission = QtWidgets.QPushButton('Full Mission', self)
        time_mission.setFixedWidth(self.default_button_width)
        time_mission.setFont(self.default_font)
        self.grid.addWidget(time_mission, self.gui_row, 1, 1, 5)
        time_mission.clicked.connect(self.event_time_mission)
        self.gui_row += 1

        search_start = QtWidgets.QPushButton('Search', self)
        search_start.setFixedWidth(self.default_button_width)
        search_start.setFont(self.default_font)
        self.grid.addWidget(search_start, self.gui_row, 0, 1, 6)
        search_start.clicked.connect(self.event_search)

        self.gui_row -= 4

    def event_time_recent(self):
        end_time = datetime.utcnow()
        start_time = end_time - timedelta(weeks=3)
        fmt = '%Y-%m-%d' #'%d-%b-%Y'
        start_time_s = start_time.strftime(fmt+' 00:00')
        end_time_s = end_time.strftime(fmt+' 23:59')
        self.start_time.setText(start_time_s)
        self.end_time.setText(end_time_s)

    def event_time_mission(self):
        fmt = '%Y-%m-%d' #'%d-%b-%Y'
        start_time_s = '2006-10-21 00:00'
        end_time_s = datetime.utcnow().strftime(fmt+' 23:59')
        self.start_time.setText(start_time_s)
        self.end_time.setText(end_time_s)

    def select_primary(self):
        """Set primary search criterion"""
        title = QtWidgets.QLabel(self)
        title.setText('Primary Search Criteria')
        title.setAlignment(QtCore.Qt.AlignBottom)
        title.setFont(self.default_font)
        self.primary_box = QtWidgets.QComboBox()
        self.primary_box.addItems([item for item in self.criteria.keys()])
        self.primary_box.setFixedWidth(self.default_button_width) # or 150
        self.primary_box.setFont(self.default_font)
        self.primary_text = QtWidgets.QLineEdit(self)
        self.primary_text.setFont(self.default_font)

        self.grid.addWidget(title, self.gui_row, 2, 1, 2)
        self.gui_row += 1
        self.grid.addWidget(self.primary_box, self.gui_row, 2)
        self.grid.addWidget(self.primary_text, self.gui_row, 3, 1, 3)

        # Advance gui row as needed
        self.gui_row += 1

    def set_filters(self):
        """Set search result filters"""
        title = QtWidgets.QLabel(self)
        title.setText('Result Filters')
        title.setAlignment(QtCore.Qt.AlignBottom)
        title.setFont(self.default_font)

        rast_type_title = QtWidgets.QLabel(self)
        rast_type_title.setText('Raster Type (#)')
        rast_type_title.setFont(self.default_font)
        self.rast_type_box = QtWidgets.QComboBox()
        self.rast_type_box.addItems([item for item in self.rast_types.keys()])
        self.rast_type_box.setFixedWidth(self.default_button_width) # or 150
        self.rast_type_box.setFont(self.default_font)

        slit_slot_title = QtWidgets.QLabel(self)
        slit_slot_title.setText('Slit/Slot (slit_index)')
        slit_slot_title.setFont(self.default_font)
        self.slit_slot_box = QtWidgets.QComboBox()
        self.slit_slot_box.addItems([item for item in self.slit_slot.keys()])
        self.slit_slot_box.setFixedWidth(self.default_button_width) # or 150
        self.slit_slot_box.setFont(self.default_font)

        wave_title = QtWidgets.QLabel(self)
        wave_title.setText(u'Wavelength(s) [\u212B]')
        wave_title.setFont(self.default_font)
        self.wave_text = QtWidgets.QLineEdit(self)
        self.wave_text.setFixedWidth(self.default_button_width) # or 150
        self.wave_text.setFont(self.default_font)

        apply_filter = QtWidgets.QPushButton('Apply Filters', self)
        apply_filter.setFixedWidth(self.default_button_width)
        apply_filter.setFont(self.default_font)
        apply_filter.clicked.connect(self.event_apply_filter)

        clear_filter = QtWidgets.QPushButton('Clear Filters', self)
        clear_filter.setFixedWidth(self.default_button_width)
        clear_filter.setFont(self.default_font)
        clear_filter.clicked.connect(self.event_clear_filter)

        self.grid.addWidget(title, self.gui_row, 2, 1, 2)
        self.gui_row += 1
        self.grid.addWidget(clear_filter, self.gui_row, 2)
        self.grid.addWidget(rast_type_title, self.gui_row, 3)
        self.grid.addWidget(slit_slot_title, self.gui_row, 4)
        self.grid.addWidget(wave_title, self.gui_row, 5)
        self.gui_row += 1
        self.grid.addWidget(apply_filter, self.gui_row, 2)
        self.grid.addWidget(self.rast_type_box, self.gui_row, 3)
        self.grid.addWidget(self.slit_slot_box, self.gui_row, 4)
        self.grid.addWidget(self.wave_text, self.gui_row, 5)

        # Advance gui row as needed
        self.gui_row += 1

    def event_search(self):
        """Validate and process search request."""
        self.tabs.setCurrentIndex(0) #switch to details tab
        self.search_info.setText('Found ?? search results')
        self.filter_info.setText('Showing ?? filter matches')
        self.info_detail.clear()
        if self.db_loaded == False:
            self.info_detail.append('No EIS As-Run Catalog found!\n\n'
                                    +'Please use the "Update Database" '
                                    +'button above.')
            return
        else:
            self.info_detail.append('Searching database. Please wait...')
        self.table_m.clearContents()
        self.table_m.setRowCount(1)
        QtWidgets.QApplication.processEvents() # update gui while user waits

        # Get dates and user input text
        start_time = str(self.start_time.text())
        end_time = str(self.end_time.text())
        primary_key = str(self.primary_box.currentText())
        primary_value = str(self.primary_text.text())

        if self.criteria[primary_key] == 'date_obs':
            self.d.get_by_date(start_time, end_time) # searches eis_experiment
        elif primary_key.lower().startswith('trigger'):
            # EIS triggered studies (1=XRT flare, 3=EIS flare, 4=EIS BP)
            search_kwargs = {'date':[start_time, end_time]}
            search_kwargs['tl_id'] = ['1', '3', '4']
            self.d.search(**search_kwargs, quiet=True)
        else:
            search_kwargs = {'date':[start_time, end_time]}
            search_kwargs[self.criteria[primary_key]] = primary_value
            self.d.search(**search_kwargs, quiet=True)

        self.selected_file = None
        if len(self.d.eis_str) > 0:
            info = []
            i = 0
            for row in self.d.eis_str:
                info.append([row.date_obs, row.study_id, row.stud_acr,
                             row.obstitle, row.xcen, row.ycen,
                             row.filename, row.tl_id, row.rastertype,
                             row.slit_index, row.wavemin, row.wavemax])
            info.sort(key=lambda x: x[6]) # sort by filename
            self.count_results = len(info)
            self.search_info.setText('Found '+str(len(info))+' search results')
            self.info_detail.append('Search complete!')
            self.info_detail.append('\nSelect any item in a row to see more'
                                   +' information')
            self.table_info = info
            self.mk_table(info)
        else:
            self.file_list = []
            self.table_info = [(None, None, None, None, None, None,
                                None, None, None, None, None, None)]
            self.count_results = 0
            self.search_info.setText('Found 0 search results')
            self.filter_info.setText('Showing 0 filter matches')
            self.info_detail.clear()
            self.info_detail.append('No entries found')

    def catalog_table(self):
        """Table with summary of search results"""
        self.search_info = QtWidgets.QLabel(self)
        self.search_info.setText('Found ?? search results')
        self.search_info.setFont(self.default_font)
        self.grid.addWidget(self.search_info, self.gui_row, 0, 1, 2)

        self.filter_info = QtWidgets.QLabel(self)
        self.filter_info.setText('Showing ?? filter matches')
        self.filter_info.setFont(self.default_font)
        self.grid.addWidget(self.filter_info, self.gui_row, 2, 1, 2)
        self.gui_row += 1

        headers = ['Date Observed', 'Study ID', 'Study Acronym',
                   'Obs. Title', 'Xcen', 'Ycen']
        widths = [160, 60, 160, 300, 60, 60] #[180, 80, 180, 350, 80, 80]

        self.table_m = QtWidgets.QTableWidget(self)
        self.table_m.verticalHeader().setVisible(False)
        self.table_m.setRowCount(1)
        self.table_m.setColumnCount(len(headers))
        self.table_m.setItem(0, 0, QtWidgets.QTableWidgetItem(' '))
        self.table_m.setHorizontalHeaderLabels(headers)
        for i in range(len(headers)):
            self.table_m.setColumnWidth(i, widths[i])
        self.grid.addWidget(self.table_m, self.gui_row, 0, 1, 6)
        self.grid.setRowStretch(self.gui_row, 1)
        self.gui_row += 1

    def mk_table(self, info):
        """Add entries to the results table."""
        len_info = len(info)
        self.file_list = []
        self.table_m.clearContents()
        self.table_m.setRowCount(1)
        r_type = self.rast_types[self.rast_type_box.currentText()]
        s_index = self.slit_slot[self.slit_slot_box.currentText()]
        wave_list = list(str(self.wave_text.text()).strip().split(','))
        self.count_filtered = 0
        for row in range(len_info):
            # Apply result filters (raster type, slit index, wavelengths)
            if r_type is not None and int(info[row][8]) != r_type:
                continue
            elif s_index is not None and int(info[row][9]) not in s_index:
                continue
            elif any(wave_list):
                missing_wave = False
                for w in range(len(wave_list)):
                    try:
                        wvl = float(wave_list[w])
                        # Note: row[10] == wavemin array, row[11] == wavemax array
                        wave_check = (wvl - info[row][10])*(info[row][11] - wvl)
                        if wave_check.max() < 0:
                            missing_wave = True
                            break
                    except:
                        # Might be good to print a warning about invalid inputs
                        pass
                if missing_wave:
                    continue

            # If row passes all filters, extend the table and append data
            new_row_ind = self.count_filtered
            self.count_filtered += 1
            self.table_m.setRowCount(self.count_filtered)

            # Date and start time
            item = QtWidgets.QTableWidgetItem(info[row][0])
            self.table_m.setItem(new_row_ind, 0, item)
            # Study ID
            item = QtWidgets.QTableWidgetItem(str(info[row][1]))
            self.table_m.setItem(new_row_ind, 1, item)
            # Study acronym
            item = QtWidgets.QTableWidgetItem(info[row][2])
            self.table_m.setItem(new_row_ind, 2, item)
            # Description
            item = QtWidgets.QTableWidgetItem(info[row][3])
            self.table_m.setItem(new_row_ind, 3, item)
            # Xcen and Ycen
            fstring = '{:0.1f}'.format(info[row][4])
            item = QtWidgets.QTableWidgetItem(fstring)
            self.table_m.setItem(new_row_ind, 4, item)
            fstring = '{:0.1f}'.format(info[row][5])
            item = QtWidgets.QTableWidgetItem(fstring)
            self.table_m.setItem(new_row_ind, 5, item)
            self.file_list.append(info[row][6])

        # Update filter count label
        if r_type is None and s_index is None and not any(wave_list):
            self.filter_info.setText('Showing all results (no filter applied)')
        else:
            self.filter_info.setText('Showing '+str(self.count_filtered)
                                    +' filter matches')

        # Any cells highlighted?
        self.table_m.currentCellChanged.connect(self.get_details)
        self.tabs.currentChanged.connect(self.event_update_context_image)

    def get_details(self, row, column):
        """Provide details on selected cell."""
        self.info_detail.clear()
        # Checking for a valid filename prevents crashes when a new search is
        # made and the results list is empty
        try:
            row_filename = str(self.file_list[row])
        except:
            row_filename = None
        if row_filename:
            info = self.fill_info(row_filename)
            for line in info:
                self.info_detail.append(line)
            self.info_detail.verticalScrollBar().\
                setValue(self.info_detail.verticalScrollBar().minimum())

            # Update the context image
            self.event_update_context_image(0)

    def fill_info(self, file):
        """Retrieve useful info for a selected file."""
        info = []
        if len(self.d.eis_str) != 0:
            row, = [x for x in self.d.eis_str if x.filename == str(file)]
            self.selected_file = row.filename
            info.append(f"{'filename':<20} {row.filename}")
            info.append(f"{'date_obs':<20} {row.date_obs}")
            info.append(f"{'date_end':<20} {row.date_end}")
            info.append(f"{'xcen':<20} {row.xcen}")
            info.append(f"{'ycen':<20} {row.ycen}")
            info.append(f"{'fovx':<20} {row.fovx}")
            info.append(f"{'fovy':<20} {row.fovy}")
            info.append(f"{'tl_id':<20} {row.tl_id}")
            info.append(f"{'study_id':<20} {row.study_id}")
            info.append(f"{'stud_acr':<20} {row.stud_acr}")
            info.append(f"{'rast_acr':<20} {row.rast_acr}")
            info.append(f"{'rast_id':<20} {row.rast_id}")
            info.append(f"{'jop_id':<20} {row.jop_id}")
            info.append(f"{'obstitle':<20} {row.obstitle}")
            info.append(f"{'obs_dec':<20} {row.obs_dec}")
            info.append(f"{'sci_obj':<20} {row.sci_obj}")
            info.append(f"{'target':<20} {row.target}")
            info.append(f"{'rastertype':<20} {self.rtype_dict[int(row.rastertype)]}")
            info.append(f"{'slit_index':<20} {self.sindex_dict[int(row.slit_index)]}")
            info.append(f"{'scan_fm_nsteps':<20} {row.scan_fm_nsteps}")
            info.append(f"{'scan_fm_stepsize':<20} {row.scan_fm_stepsize}")
            info.append(f"{'nexp':<20} {row.nexp}")
            info.append(f"{'exptime':<20} {row.exptime}")

            info.append(f"\n\n{'----- Line List -----':^55}")
            info.append(f"{'window':<8} {'title':<20} "
                       +f"{'wavemin':<9} {'wavemax':<9} {'width':<5}")
            for i in range(0,len(row.ll_title)):
                info.append(f"{i:<8} {row.ll_title[i]:<20} "
                           +f"{row.wavemin[i]:<9.2f} {row.wavemax[i]:<9.2f} "
                           +f"{row.width[i]:<5}")
            info.append("\n")

            # Generate info for display in context tab
            self.selected_info = []
            self.selected_info.append(f"{'filename':<20} {row.filename}")
            self.selected_info.append(f"{'date_obs, date_end':<20} {row.date_obs}   to   {row.date_end}")
            self.selected_info.append(f"{'xcen, ycen':<20} {row.xcen:0.2f}, {row.ycen:0.2f}")
            self.selected_info.append(f"{'fovx, fovy':<20} {row.fovx:0.2f}, {row.fovy:0.2f}")
            self.selected_info.append(f"{'study_id, stud_acr':<20} {row.study_id}, {row.stud_acr}")
            self.selected_info.append(f"{'obstitle':<20} {row.obstitle}")
        return info

    @QtCore.pyqtSlot(int)
    def get_image(self):
        url = self.context_url
        self.start_request(url)

    def start_request(self, url):
        request = QtNetwork.QNetworkRequest(QtCore.QUrl(url))
        self.manager.get(request)

    @QtCore.pyqtSlot(QtNetwork.QNetworkReply)
    def on_finished(self, reply):
        target = reply.attribute(QtNetwork.QNetworkRequest.RedirectionTargetAttribute)
        if reply.error():
            print("error: {}".format(reply.errorString()))
            return
        elif target:
            newUrl = reply.url().resolved(target)
            self.start_request(newUrl)
            return
        pixmap = QtGui.QPixmap()
        pixmap.loadFromData(reply.readAll())
        pixmap = pixmap.scaled(self.context_imgNX, self.context_imgNY)
        self.context_img.setPixmap(pixmap)

    def event_update_context_image(self, tab_index):
        """Download context image into memory and update the image tab"""
        clean_filename = self.selected_file.replace('.gz', '')
        remote_dir = get_remote_image_dir(clean_filename)
        context_img_name = 'XRT_'+clean_filename+'.gif'
        self.context_url = remote_dir+context_img_name

        if self.tabs.currentIndex() == 1:
            try:
                self.get_image()
            except:
                print('   ERROR: context images or server are unavailable.')

            self.info_context.clear()
            for line in self.selected_info:
                self.info_context.append(line)
        else:
            self.event_clear_context_image()

    def event_clear_context_image(self):
        self.info_context.clear()
        self.info_context.append('Select any item in a row to see a context'
                                +' image')
        buff = np.zeros((self.context_imgNX, self.context_imgNX, 3), dtype=np.int16)
        image = QtGui.QImage(buff, self.context_imgNX, self.context_imgNY, QtGui.QImage.Format_ARGB32)
        self.context_img.setPixmap(QtGui.QPixmap(image))

    def details(self):
        """Display detailed cat info."""
        # Initialize tab panel and add main tabs
        self.tabs = QtWidgets.QTabWidget()
        self.detail_tab = QtWidgets.QWidget()
        self.image_tab = QtWidgets.QWidget()
        self.help_tab = QtWidgets.QWidget()
        self.tabs.addTab(self.detail_tab,"Details")
        self.tabs.addTab(self.image_tab,"Images")
        self.tabs.addTab(self.help_tab,"Help")

        # Create details tab
        self.detail_tab.grid = QtWidgets.QGridLayout()
        self.info_detail = QtWidgets.QTextEdit()
        self.info_detail.setFont(self.info_detail_font)
        self.info_detail.setReadOnly(True)
        self.detail_tab.grid.addWidget(self.info_detail)
        self.detail_tab.setLayout(self.detail_tab.grid)

        # Create the image tab and initialize SSL context for downloading
        self.image_tab.grid = QtWidgets.QGridLayout()
        self.info_context = QtWidgets.QTextEdit()
        self.info_context.setFont(self.info_detail_font)
        self.info_context.setReadOnly(True)
        self.image_tab.grid.addWidget(self.info_context, 0, 0)
        self.context_img = QtWidgets.QLabel()
        buff = np.zeros((self.context_imgNX, self.context_imgNX, 3), dtype=np.int16)
        image = QtGui.QImage(buff, self.context_imgNX, self.context_imgNY,
                             QtGui.QImage.Format_ARGB32)
        self.context_img.setPixmap(QtGui.QPixmap(image))
        self.image_tab.grid.addWidget(self.context_img, 1, 0, 4, 1)
        self.image_tab.setLayout(self.image_tab.grid)
        # self.ssl_context = ssl.SSLContext(ssl.PROTOCOL_TLS)
        # self.ssl_context.load_verify_locations(certifi.where())

        # Create help tab
        self.help_tab.grid = QtWidgets.QGridLayout()
        self.info_help = QtWidgets.QTextEdit()
        self.info_help.setFont(self.info_detail_font)
        self.info_help.setReadOnly(True)
        self.help_tab.grid.addWidget(self.info_help)
        self.help_tab.setLayout(self.help_tab.grid)

        # Add tabs to main window
        self.tabs.setStyleSheet('QTabBar{font-size: 11pt; font-family: Courier New;}')
        self.grid.addWidget(self.tabs, 1, 6, self.gui_row + 1, 3)
        self.tabs.setCurrentIndex(2) # switch to help tab

    def event_apply_filter(self):
        if self.count_results > 0:
            self.mk_table(self.table_info)

    def event_clear_filter(self):
        self.rast_type_box.setCurrentIndex(0)
        self.slit_slot_box.setCurrentIndex(0)
        self.wave_text.clear()
        if self.count_results > 0:
            self.mk_table(self.table_info)

    def save_options(self):
        """Controls for saving files."""

        data_source_title = QtWidgets.QLabel(self)
        data_source_title.setText('Data Source')
        data_source_title.setFont(self.default_font)
        self.grid.addWidget(data_source_title, self.gui_row, 0)

        set_save_dir = QtWidgets.QPushButton('Change Save Dir', self)
        set_save_dir.setFixedWidth(self.default_button_width)
        set_save_dir.setFont(self.default_font)
        self.grid.addWidget(set_save_dir, self.gui_row, 1)
        set_save_dir.clicked.connect(self.event_set_save_dir)

        self.topdir_box = QtWidgets.QLineEdit(self)
        self.topdir_box.resize(4*self.default_button_width, self.frameGeometry().height())
        self.topdir_box.setText(self.default_topdir)
        self.topdir_box.setFont(self.default_font)
        self.grid.addWidget(self.topdir_box, self.gui_row, 2, 1, 4)

        self.gui_row += 1

        self.data_source_box = QtWidgets.QComboBox()
        self.data_source_box.addItems(['NRL (main)', 'NASA (mirror)', 'MSSL (mirror)'])
        self.data_source_box.setFixedWidth(self.default_button_width) # or 150
        self.data_source_box.setFont(self.default_font)
        self.grid.addWidget(self.data_source_box, self.gui_row, 0)

        download_selected = QtWidgets.QPushButton('Download Selected', self)
        download_selected.setFixedWidth(self.default_button_width)
        download_selected.setFont(self.default_font)
        self.grid.addWidget(download_selected, self.gui_row, 1)
        download_selected.clicked.connect(self.event_download_selected)

        download_list = QtWidgets.QPushButton('Download All', self)
        download_list.setFixedWidth(self.default_button_width)
        download_list.setFont(self.default_font)
        self.grid.addWidget(download_list, self.gui_row, 2)
        download_list.clicked.connect(self.event_download_file_list)

        self.radio = QtWidgets.QRadioButton("Use Date Tree")
        self.radio.setFixedWidth(self.default_button_width)
        self.radio.setFont(self.default_font)
        self.grid.addWidget(self.radio, self.gui_row, 3)

        self.save_list = QtWidgets.QPushButton('Save File List', self)
        self.save_list.setFixedWidth(self.default_button_width)
        self.save_list.setFont(self.default_font)
        self.grid.addWidget(self.save_list, self.gui_row, 4)
        self.save_list.clicked.connect(self.event_save_file_list)

        self.filename_box = QtWidgets.QLineEdit(self)
        self.filename_box.setFixedWidth(self.default_button_width)
        self.filename_box.setText(self.default_filename)
        self.filename_box.setFont(self.default_font)
        self.grid.addWidget(self.filename_box, self.gui_row, 5)

    def event_set_save_dir(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        options |= QtWidgets.QFileDialog.ShowDirsOnly
        new_savedir = QtWidgets.QFileDialog.getExistingDirectory(self,
                                                'Select a directory',
                                                options=options)
        if os.path.isdir(new_savedir):
            self.topdir_box.setText(new_savedir)

    def event_download_selected(self):
        if self.selected_file is not None:
            self.tabs.setCurrentIndex(0) #switch to details tab
            data_source = self.data_source_box.currentText()
            datetree = self.radio.isChecked()
            topdir = self.topdir_box.text()
            self.info_detail.clear()
            info = (f'Downloading {self.selected_file}\n'
                    f'   Data Source: {data_source}\n'
                    f'   Save dir: {topdir}\n\n'
                    f'Please wait...')
            self.info_detail.append(info)
            QtWidgets.QApplication.processEvents() # update gui while user waits
            o = download_hdf5_data(filename=self.selected_file, datetree=datetree,
                                   source=data_source, local_top=topdir, overwrite=True)
            self.info_detail.append('\nComplete')

    def event_download_file_list(self):
        if self.file_list is not None:
            self.tabs.setCurrentIndex(0) #switch to details tab
            data_source = self.data_source_box.currentText()
            datetree = self.radio.isChecked()
            topdir = self.topdir_box.text()
#            sys.stdout = OutLog(self.info_detail, sys.stdout)
#            sys.stderr = OutLog(self.info_detail, sys.stderr)
            self.info_detail.clear()
            info = (f'Downloading all files listed above\n'
                    f'   Data Source: {data_source}\n'
                    f'   Save dir: {topdir}\n\n'
                    f'Please wait (see console for download progress)...')
            self.info_detail.append(info)
            QtWidgets.QApplication.processEvents() # update gui while user waits
            o = download_hdf5_data(filename=self.file_list, datetree=datetree,
                                   source=data_source, local_top=topdir, overwrite=True)
            self.info_detail.append('\nComplete')

    def event_save_file_list(self):
        """Save a list of the displayed files, one per line."""
        if self.file_list is not None:
            topdir = self.topdir_box.text()
            filename = self.filename_box.text()
            if filename != '':
                with open(os.path.join(topdir, filename), 'w') as fp:
                    for item in self.file_list:
                        fp.write("{}\n".format(item))
                        print(item)
                print(f' + wrote {filename}')

# this class handles the redirection of stdout, not implemented yet.
class OutLog:
    def __init__(self, edit, out=None):
        self.edit = edit
        self.out = None

    def write(self, m):
        self.edit.moveCursor(QtGui.QTextCursor.End)
        self.edit.insertPlainText(m)

        if self.out:
            self.out.write(m)


#-#-#-#-# MAIN #-#-#-#-#
def eis_catalog():
    if len(sys.argv) == 2:
        # option 1 - user inputs the file
        db_file = sys.argv[1]
    else:
        # option 2a - use file from SSW (fully configured environments)
        ssw_dir = os.environ.get('SSW')
        if ssw_dir is None:
            # option 2b - use file from SSW (search for SSW installation)
            print('NOTICE: "SSW" environment variable not found.')
            print('   Searching common installation directories...')
            for test_dir in ['/usr/local/ssw', os.path.expanduser('~')+'/ssw',
                             'C:\\ssw', 'D:\\ssw']:
                if os.path.isdir(test_dir):
                    ssw_dir = test_dir
                    print('SSW found in '+ssw_dir)
                    break
        if ssw_dir is not None:
            db_file = os.path.join(ssw_dir, 'hinode', 'eis', 'database',
                                   'catalog', 'eis_cat.sqlite')
        else:
            # option 3 - use local database stored within eispac
            print('WARNING: No SSW installation found.')
            print('   Attempting to use local eis database.')
            import eispac
            module_dir = os.path.dirname(eispac.download.__file__)
            db_file = os.path.join(module_dir, 'eis_cat.sqlite')

    app = QtWidgets.QApplication([])
    topthing = Top(db_file)

    sys.exit(app.exec_())

if __name__ == '__main__':
    eis_catalog()

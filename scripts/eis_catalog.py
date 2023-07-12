#!/usr/bin/env python
"""EIS asrun catalog search PyQt4/5 gui.

A python gui to search the EIS asrun catalog. The intention is to
be somewhat like the IDL routine eis_cat.pro, but with fewer
features, at least initially. Very much a work in progress.

(2017-Apr-12) First working version added to git.
(2020-Dec-16) Now tries to find SSW if the 'SSW' environment variable is missing
(2020-Dec-18) If no database is found, ask the user if they want to download it.
(2021-Oct-22) Major update adding in more search criteria and faster queries.
"""
__all__ = ['eis_catalog']

import sys
import os, time
import re
import sqlite3
import pkgutil
from datetime import datetime, timedelta
from eispac.download import eis_obs_struct
from eispac.download.download_hdf5_data import download_hdf5_data
from eispac.download.download_db import download_db

if pkgutil.find_loader('PyQt4'):
    from PyQt4 import QtCore
    from PyQt4 import QtGui as QtWidgets
    have_Qt4 = True
    have_Qt5 = False
elif pkgutil.find_loader('PyQt5'):
    from PyQt5 import QtCore, QtWidgets, QtGui
    have_Qt4 = False
    have_Qt5 = True
else:
    print('Either PyQt4 or PyQt5 must be installed.')
    sys.exit()

class Top(QtWidgets.QWidget):

    def __init__(self, dbfile, parent=None):
        super(Top, self).__init__(parent)
        self.file_list = None
        self.selected_file = None
        self.default_filename = 'eis_filelist.txt'
        self.default_start_time = '2018-05-29 00:00' # '29-May-2018 00:00'
        self.default_end_time = '2018-05-29 23:59' # '29-May-2018 23:59'
        self.default_button_width = 165 #130
        self.default_topdir = os.path.join(os.getcwd(), 'data_eis')
        self.dbfile = dbfile
        self.db_loaded = False

        # Font settings
        if have_Qt4:
            self.default_font = QtWidgets.QFont()
            self.small_font = QtWidgets.QFont()
            self.info_detail_font = QtWidgets.QFont("Courier New", 9)
        else:
            self.default_font = QtGui.QFont()
            self.small_font = QtGui.QFont()
            self.info_detail_font = QtGui.QFont("Courier New", 9)
        self.default_font.setPointSize(11)
        self.small_font.setPointSize(9)

        # Dict of search criteria to include as options in the drop-down list
        # Note: the keys:value pairs give the mapping of gui_label:sqlite_col
        self.criteria = {'Date Only':'date_obs', 'Study ID':'study_id',
                         'Study Acronym':'stud_acr', 'HOP ID':'jop_id',
                         'Target':'target', 'Sci. Obj.':'sci_obj',
                         'Obs. Title':'obstitle'}
        # Filter drop-down lists. given as gui_label:filter_value pairs
        self.rast_types = {'Any':None, 'Scan (0)':0, 'Sit-and-Stare (1)':1}
        self.slit_slot = {'Any':None, '1" slit (0)':0, '2" slit (2)':2,
                          '40" slot (3)':3, '266" slot (1)':1}

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

        self.init_ui()

    def init_ui(self):
        """Manage everything."""
        self.grid = QtWidgets.QGridLayout(self)
        self.gui_row = 0
        # self.setStyleSheet("QLabel{font-size: 11pt;}")

        # Quit, Help, Update DB buttons (with filename & timestamp)
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
        self.setGeometry(50, 100, 1800, 800)
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

        self.help = QtWidgets.QPushButton('Help', self)
        self.help.setFixedWidth(self.default_button_width)
        self.help.setFont(self.default_font)
        self.help.clicked.connect(self.event_help)

        self.download_db = QtWidgets.QPushButton('Update Database', self)
        self.download_db.setFixedWidth(self.default_button_width)
        self.download_db.setFont(self.default_font)
        self.download_db.clicked.connect(self.event_download_db)

        self.db_source_box = QtWidgets.QComboBox()
        self.db_source_box.addItems(['Hesperia (source)', 'NRL (mirror)'])
        self.db_source_box.setFixedWidth(self.default_button_width) # or 150
        self.db_source_box.setFont(self.default_font)

        self.db_info = QtWidgets.QLabel(self)
        self.db_info.setFixedWidth(3*self.default_button_width)
        self.db_info.setFont(self.small_font)
        if os.path.isfile(self.dbfile):
            self.update_db_file_label()
        else:
            self.db_info.setText('Unable to locate DB: ' + self.dbfile)

        self.grid.addWidget(self.quit, self.gui_row, 0)
        self.grid.addWidget(self.help, self.gui_row, 1)
        self.grid.addWidget(self.download_db, self.gui_row, 2)
        self.grid.addWidget(self.db_source_box, self.gui_row, 3)
        self.grid.addWidget(self.db_info, self.gui_row, 4, 1, 3)

        self.gui_row += 1

    def update_db_file_label(self):
         t = 'DB file: '+os.path.abspath(self.dbfile)+'\nDownload date: '+ \
                time.ctime(os.path.getmtime(self.dbfile)).lstrip().rstrip()
         self.db_info.setText(t)

    def event_download_db(self):
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
        self.info_detail.clear()
        help_text = """EIS As-Run Catalog Search Tool

* Please use ISO format dates (e.g., YYYY-MM-DD HH:MM). If only
  the start date is provided, the end date will be assumed to be
  24 hours later. WARNING: "Date Only" searches over the entire
  mission can take a VERY long time, please be patient.

* Primary search criteria descriptions:

  Study ID
      ID number for specific observation plan (line list, raster
      steps, etc.). Studies may repeated throughout the mission.

  Study Acronym
      Short text label for the study. You only need to input the
      first few characters (case is ignored).

  HOP ID
      Hinode Operation Plan ID. Assigned to observations that were
      coordinated with other telescopes or spacecraft missions.
      Also known as "JOP ID" (Joint Obs. Program ID)

  Target
      Main observation target (e.g., Active Region, Quiet Sun).

  Sci. Obj.
      Target phenomena (e.g., AR, QS, BP, LMB)

  Obs. Title
      Observation title. Just a word or two is enough, results will
      be shown for all titles containing the input text.

  Please Note: If "Target" and "Sci. Obj." are not defined by the
  study author, default values "Quiet Sun" & "QS" are assigned,
  regardless of the actual observation target.

* If the "Use Date Tree" box is checked, files will be downloaded
  into subdirectories organized by date (../YYYY/MM/DD/)
"""
        self.info_detail.append(help_text)
        self.info_detail.verticalScrollBar().setValue(
                                self.info_detail.verticalScrollBar().minimum())

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
            self.d.query_main(date=[start_time, end_time])
        else:
            search_kwargs = {'date':[start_time, end_time]}
            search_kwargs[self.criteria[primary_key]] = primary_value
            self.d.query_main(**search_kwargs)

        self.selected_file = None
        if len(self.d.eis_str) > 0:
            # print('nrows: ', len(self.d.eis_str))
            # print(eis_obs_struct.members(self.d.eis_str[0]))
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
        # title = QtWidgets.QLabel(self)
        # title.setText('Catalog Search Results')
        # title.setFont(self.default_font)
        # self.grid.addWidget(title, self.gui_row, 0, 1, 2)

        self.search_info = QtWidgets.QLabel(self)
        self.search_info.setText('Found ?? search results')
        self.search_info.setFont(self.default_font)
        self.grid.addWidget(self.search_info, self.gui_row, 0, 1, 2)

        self.filter_info = QtWidgets.QLabel(self)
        self.filter_info.setText('Showing ?? filter matches')
        self.filter_info.setFont(self.default_font)
        self.grid.addWidget(self.filter_info, self.gui_row, 2, 1, 2)
        self.gui_row += 1

        # headers = ['Date Observed', 'Study ID', 'Study Acronym',
        #            'Obs. Title', 'Xcen', 'Ycen', 'Filename']
        # widths = [180, 80, 180, 350, 80, 80, 210]

        headers = ['Date Observed', 'Study ID', 'Study Acronym',
                   'Obs. Title', 'Xcen', 'Ycen']
        widths = [180, 80, 180, 350, 80, 80]

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
        # self.table_m.setRowCount(len_info)
        r_type = self.rast_types[self.rast_type_box.currentText()]
        s_index = self.slit_slot[self.slit_slot_box.currentText()]
        wave_list = list(str(self.wave_text.text()).strip().split(','))
        self.count_filtered = 0
        for row in range(len_info):
            # Apply result filters (raster type, slit index, wavelengths)
            if r_type is not None and int(info[row][8]) != r_type:
                continue
            elif s_index is not None and int(info[row][9]) != s_index:
                continue
            elif any(wave_list):
                missing_wave = False
                for w in range(len(wave_list)):
                    try:
                        wvl = float(wave_list[w])
                        # row[10] == wavemin array, row[11] == wavemax array
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
            # # Filename
            # item = QtWidgets.QTableWidgetItem(info[row][6])
            # self.table_m.setItem(row, 6, item)
            self.file_list.append(info[row][6])

        # Update filter count label
        if r_type is None and s_index is None and not any(wave_list):
            self.filter_info.setText('Showing all results (no filter applied)')
        else:
            self.filter_info.setText('Showing '+str(self.count_filtered)
                                    +' filter matches')

        # Any cells highlighted?
        self.table_m.cellClicked.connect(self.get_details)

    def get_details(self, row, column):
        """Provide details on selected cell."""
        self.info_detail.clear()
        # Checking for a valid filename prevents crashes when a new search is
        # made and the results list is empty
        try:
            # row_filename = str(self.table_info[row][6])
            row_filename = str(self.file_list[row])
        except:
            row_filename = None
        if row_filename:
            info = self.fill_info(row_filename)
            for line in info:
                self.info_detail.append(line)
            self.info_detail.verticalScrollBar().\
                setValue(self.info_detail.verticalScrollBar().minimum())

    def fill_info(self, file):
        """Retrieve useful info for a selected file."""
        info = []
        if len(self.d.eis_str) != 0:
            row, = [x for x in self.d.eis_str if x.filename == str(file)]
            self.selected_file = row.filename
            info.append("{0:<20} {1}".format('filename', row.filename))
            info.append("{0:<20} {1}".format('date_obs', row.date_obs))
            info.append("{0:<20} {1}".format('date_end', row.date_end))
            info.append("{0:<20} {1:0.1f}".format('xcen', row.xcen))
            info.append("{0:<20} {1:0.1f}".format('ycen', row.ycen))
            info.append("{0:<20} {1:0.1f}".format('fovx', row.fovx))
            info.append("{0:<20} {1:0.1f}".format('fovy', row.fovy))
            info.append("{0:<20} {1}".format('tl_id', row.tl_id))
            info.append("{0:<20} {1}".format('study_id', row.study_id))
            info.append("{0:<20} {1}".format('stud_acr', row.stud_acr))
            info.append("{0:<20} {1}".format('rast_acr', row.rast_acr))
            info.append("{0:<20} {1}".format('rast_id', row.rast_id))
            info.append("{0:<20} {1}".format('jop_id', row.jop_id))
            info.append("{0:<20} {1}".format('obstitle', row.obstitle))
            info.append("{0:<20} {1}".format('obs_dec', row.obs_dec))
            info.append("{0:<20} {1}".format('sci_obj', row.sci_obj))
            info.append("{0:<20} {1}".format('target', row.target))
            info.append("{0:<20} {1}".format('rastertype', row.rastertype))
            info.append("{0:<20} {1}".format('slit_index', row.slit_index))
            info.append("{0:<20} {1}".format('scan_fm_nsteps', row.scan_fm_nsteps))
            info.append("{0:<20} {1}".format('scan_fm_stepsize',
                                             row.scan_fm_stepsize))
            info.append("{0:<20} {1}".format('nexp', row.nexp))
            info.append("{0:<20} {1}".format('exptime', row.exptime))

            info.append(f"\n\n{'----- Line List -----':^55}")
            info.append(f"{'window':<8} {'title':<20} "
                       +f"{'wavemin':<9} {'wavemax':<9} {'width':<5}")
            for i in range(0,len(row.ll_title)):
                info.append(f"{i:<8} {row.ll_title[i]:<20} "
                           +f"{row.wavemin[i]:<9.2f} {row.wavemax[i]:<9.2f} "
                           +f"{row.width[i]:<5}")
            info.append("\n")
        return info

    def details(self):
        """Display detailed cat info."""
        title = QtWidgets.QLabel(self)
        title.setText('Details')
        title.setFont(self.default_font)
        # self.grid.addWidget(title, self.gui_row, 0, 1, 6)
        # self.grid.addWidget(title, self.gui_row, 6, 1, 3)
        self.grid.addWidget(title, 1, 6, 1, 3)
        # self.gui_row += 1

        self.info_detail = QtWidgets.QTextEdit()
        self.info_detail.setFont(self.info_detail_font)
        # self.grid.addWidget(self.info_detail, self.gui_row, 0, 1, 6)
        # self.grid.addWidget(self.info_detail, self.gui_row, 6, 1, 3)
        self.grid.addWidget(self.info_detail, 2, 6, self.gui_row, 3)
        #self.grid.setRowStretch(self.gui_row, 1)
        # self.gui_row += 1
        self.info_detail.setReadOnly(True)

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
        # self.topdir_box.setFixedWidth(self.default_button_width)
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

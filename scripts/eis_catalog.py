#!/usr/bin/env python
"""EIS asrun catalog search PyQt4/5 gui.

A python gui to search the EIS asrun catalog. The intention is to
be somewhat like the IDL routine eis_cat.pro, but with fewer
features, at least initially. Very much a work in progress.

(2017-Apr-12) First working version added to git.
(2020-Dec-16) Now tries to find SSW if the 'SSW' environment variable is missing
(2020-Dec-18) If no database is found, as the user if they want to download it.

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
        self.default_filename = 'eis_filelist.txt'
        self.default_start_time = '2018-05-29 00:00' # '29-May-2018 00:00'
        self.default_end_time = '2018-05-29 23:59' # '29-May-2018 23:59'
        self.default_button_width = 130
        self.default_topdir = 'data_eis/'
        self.dbfile =  dbfile

        if os.path.isfile(self.dbfile):
            self.d = eis_obs_struct.EIS_DB(self.dbfile)
        else:
            # Ask if the user wants to download a copy of the database
            ask_db = QtWidgets.QMessageBox.question(self, 'Missing database',
                        'No EIS as-run database found\n'
                        'Would you like to download a local copy?',
                        QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
            if ask_db == QtWidgets.QMessageBox.Yes:
                self.dbfile = download_db()
                if not string(self.dbfile).endswith('.part'):
                    self.d = eis_obs_struct.EIS_DB(self.dbfile)
                else:
                    print('Failed to download EIS database! Exiting...')
                    exit()
            else:
                print('ERROR: No EIS as-run database found. Exiting...')
                exit()

        self.init_ui()

    def init_ui(self):
        """Manage everything."""
        self.grid = QtWidgets.QGridLayout(self)

        # Info at the top
        self.gui_row = 0
        self.top()

        # Search criteria
        self.search_left()
        self.search_center()
        self.search_right()

        # Catalog info
        self.main()

        # Info for one item, a possible future addition.
        self.details()

        # Bottom stuff
        self.bottom()

        # And away we go
        self.setLayout(self.grid)
        self.setGeometry(200, 200, 1240, 800)
        self.setWindowTitle('EIS As-Run Catalog Information')
        self.show()

    def event_quit(self):
        QtWidgets.QApplication.instance().quit()

    def top(self):
        """Minimal info at the top."""
        self.quit = QtWidgets.QPushButton('Quit')
        self.quit.setFixedWidth(self.default_button_width)
        self.quit.clicked.connect(self.event_quit)

        self.help = QtWidgets.QPushButton('Help', self)
        self.help.setFixedWidth(150)
        self.help.clicked.connect(self.event_help)

        self.grid.addWidget(self.quit, self.gui_row, 0)
        self.grid.addWidget(self.help, self.gui_row, 1)

        self.gui_row += 1

        self.top_title =  QtWidgets.QLabel(self)
        if os.path.isfile(self.dbfile):
            self.update_db_file_label()
        else:
            self.top_title.setText('Unable to locate file: ' + self.dbfile)
        self.grid.addWidget(self.top_title, self.gui_row, 1, 1, 5)

        self.download_db = QtWidgets.QPushButton('Get Latest DB', self)
        self.download_db.setFixedWidth(self.default_button_width)
        self.grid.addWidget(self.download_db, self.gui_row, 0)
        self.download_db.clicked.connect(self.event_download_db)

        self.gui_row += 1

    def update_db_file_label(self):
         t = 'Using file: '+os.path.abspath(self.dbfile)+'  DB Downloaded : '+ \
                time.ctime(os.path.getmtime(self.dbfile)).lstrip().rstrip()
         self.top_title.setText(t)

    def event_download_db(self):
        self.info_detail.clear()
        info = (f'Downloading eis_cat.sqlite.\n'
                f'   Remote: https://hesperia.gsfc.nasa.gov/ssw/hinode/eis'
                f'/database/catalog/eis_cat.sqlite\n'
                f'   Local: {os.path.abspath(self.dbfile)}\n\n'
                f'Please wait...')
        self.info_detail.append(info)
        QtWidgets.QApplication.processEvents() # update gui while user waits
        self.d.cur.close()
        self.d.conn.close()
        self.d = 0
        self.dbfile = download_db(os.path.dirname(self.dbfile))
        self.d = eis_obs_struct.EIS_DB(self.dbfile)
        self.update_db_file_label()
        self.info_detail.append('\nComplete')

    def event_help(self):
        """Put help info in details window."""
        self.info_detail.clear()
        info = """EIS As-Run Catalog Search Tool

The search can be on a date range (ISO standard, e.g., yyyy-mm-dd
and optionally hh:mm) only, a date range and either the Study ID
or the Study Acronym, or the Study ID or the Study Acronym alone.
If only the start date is provided, a search of that date and the
following 24 hours is assumed. Searches with no date range can
take a VERY long time, so be patient.

Highlighting any item in a row in the search result will display
more detailed information from the database about the observation
in this window.

Clicking on the Save File List button will produce a text file
with one file per line, which can be used as input to other
programs, for example, get_eis_level0.py, which will retreive the
level 0 files from ISAS.
"""
        self.info_detail.append(info)
        self.info_detail.verticalScrollBar().\
            setValue(self.info_detail.verticalScrollBar().minimum())

    def event_db_file_search(self):
        """Find a db file to use, if needed."""
        db_dialog = QtWidgets.QFileDialog(self)
        if have_Qt5:
            file_name, _ = db_dialog.getOpenFileName()
        else:
            file_name = db_dialog.getOpenFileName()
        save = self.gui_row
        self.gui_row = 0
        if file_name != '':
            self.top(file_name)
        self.gui_row = save

    def search_left(self):
        """Process search criteria, dates."""
        title = QtWidgets.QLabel(self)
        title.setText('Search Criteria. Selecting a time range ' +
                      '(ISO standard) HIGHLY recommended)')
        self.grid.addWidget(title, self.gui_row, 0, 1, 6)
        self.gui_row += 1

        start_t = QtWidgets.QLabel(self)
        start_t.setText('Start Time')
        self.start_time = QtWidgets.QLineEdit(self)
        self.start_time.setFixedWidth(self.default_button_width)
        self.start_time.setText(self.default_start_time)

        self.grid.addWidget(start_t, self.gui_row, 0)
        self.grid.addWidget(self.start_time, self.gui_row, 1)
        self.gui_row += 1

        end_t = QtWidgets.QLabel(self)
        end_t.setText('End Time')
        self.end_time = QtWidgets.QLineEdit(self)
        self.end_time.setFixedWidth(self.default_button_width)
        self.end_time.setText(self.default_end_time)
        self.grid.addWidget(end_t, self.gui_row, 0)
        self.grid.addWidget(self.end_time, self.gui_row, 1)
        self.gui_row += 1

        search_start = QtWidgets.QPushButton('Search', self)
        search_start.setFixedWidth(self.default_button_width)
        self.grid.addWidget(search_start, self.gui_row, 0, 1, 6)
        search_start.clicked.connect(self.event_search_new)

        time_set = QtWidgets.QPushButton('Last 3 Weeks', self)
        time_set.setFixedWidth(self.default_button_width)
        self.grid.addWidget(time_set, self.gui_row, 1, 1, 6)
        time_set.clicked.connect(self.event_time_set)
        self.gui_row -= 2

    def event_time_set(self):
        end_time = datetime.utcnow()
        start_time = end_time - timedelta(weeks=3)
        fmt = '%d-%b-%Y'
        start_time_s = start_time.strftime(fmt+' 00:00')
        end_time_s = end_time.strftime(fmt+' 23:59')
        self.start_time.setText(start_time_s)
        self.end_time.setText(end_time_s)

    def search_center(self):
        """Select by study ID."""
        title = QtWidgets.QLabel(self)
        title.setText('Select by Study ID')
        title.setAlignment(QtCore.Qt.AlignBottom)
        study_t = QtWidgets.QLabel(self)
        study_t.setText('Study ID')
        self.study = QtWidgets.QLineEdit(self)
        self.study.setFixedWidth(150)

        self.grid.addWidget(title, self.gui_row, 2, 1, 2)
        self.gui_row += 1
        self.grid.addWidget(study_t, self.gui_row, 2)
        self.grid.addWidget(self.study, self.gui_row, 3)
        self.gui_row += 1

        self.b1 = QtWidgets.QCheckBox('Use ID')

        self.grid.addWidget(self.b1, self.gui_row, 2, 1, 2)
        self.gui_row -= 2

    def search_right(self):
        """Select by study acronym."""
        title = QtWidgets.QLabel(self)
        title.setText('Select by Acronym (first few characters will do)')
        title.setAlignment(QtCore.Qt.AlignBottom)
        acro_t = QtWidgets.QLabel(self)
        acro_t.setText('Study Acronym')
        self.acronym = QtWidgets.QLineEdit(self)

        self.grid.addWidget(title, self.gui_row, 4, 1, 2)
        self.gui_row += 1
        self.grid.addWidget(acro_t, self.gui_row, 4)
        self.grid.addWidget(self.acronym, self.gui_row, 5)
        self.gui_row += 1

        self.b2 = QtWidgets.QCheckBox('Use Acronym')

        self.grid.addWidget(self.b2, self.gui_row, 4, 1, 2)
        self.gui_row += 1

    def event_search(self):
        """Validate and process search request."""

        # First must have dates
        start_time = str(self.start_time.text())
        end_time = str(self.end_time.text())

        self.d.get_by_date(start_time, end_time)
        if len(self.d.eis_str) > 0:
            # print('nrows: ', len(self.d.eis_str))
            # print(eis_obs_struct.members(self.d.eis_str[0]))
            info = []
            i = 0
            for row in self.d.eis_str:
                info.append([row.tl_id, row.stud_acr, row.date_obs,
                             row.obstitle, row.filename, row.xcen,
                             row.ycen, row.study_id])
        else:
            print('No entries found')

        info.sort(key=lambda x: x[2])

        if self.b1.checkState() == 2:
            s_id = str(self.study.text())
            info = [x for x in info if str(x[7]) == s_id]
        elif self.b2.checkState() == 2:
            text = str(self.acronym.text())
            r = re.compile(text, flags=re.IGNORECASE)
            info = [x for x in info if r.match(x[1])]

        self.mk_table(info)

    def event_search_new(self):
        """Validate and process search request."""
        start_time = str(self.start_time.text())
        end_time = str(self.end_time.text())

        if self.b1.checkState() == 2 and self.b2.checkState() == 2:
            self.info_detail.clear()
            self.info_detail.append("Both Study ID and Study Acronym " +
                                    "should not be checked.")
            self.info_detail.append("Using Study ID only.")

        # Roll through all the cases
        if start_time != '' and self.b1.checkState() != 2 \
                            and self.b2.checkState() != 2:
            self.d.get_by_date(start_time, end_time)
        elif start_time != '' and self.b1.checkState() == 2:
            s_id = str(self.study.text())
            self.d.get_by_study_id(s_id, date=[start_time, end_time])
        elif start_time != '' and self.b2.checkState() == 2:
            text = str(self.acronym.text())
            self.d.get_by_acronym(text, date=[start_time, end_time])
        elif start_time == '' and self.b1.checkState() == 2:
            s_id = str(self.study.text())
            self.d.get_by_study_id(s_id)
        elif start_time == '' and self.b2.checkState() == 2:
            text = str(self.acronym.text())
            self.d.get_by_acronym(text)

        if len(self.d.eis_str) > 0:
            # print('nrows: ', len(self.d.eis_str))
            # print(eis_obs_struct.members(self.d.eis_str[0]))
            info = []
            i = 0
            for row in self.d.eis_str:
                info.append([row.tl_id, row.stud_acr, row.date_obs,
                             row.obstitle, row.filename, row.xcen,
                             row.ycen, row.study_id])
            info.sort(key=lambda x: x[4]) # sort on file name
            self.mk_table(info)
        else:
            self.info_detail.clear()
            self.file_list = []
            self.info_detail.append('No entries found')

    def main(self):
        """Display main cat info."""
        title = QtWidgets.QLabel(self)
        title.setText('Catalog Search Results')
        self.grid.addWidget(title, self.gui_row, 0, 1, 2)
        self.gui_row += 1

        headers = ['Timeline ID', 'Study Acronym', 'Date Observed',
                   'Description', 'Filename', 'Xcen', 'Ycen']
        widths = [100, 180, 170, 350, 200, 80, 80]

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
        num_cols = len(info[0])
        self.file_list = []
        self.table_m.clearContents()
        self.table_m.setRowCount(len_info)
        for row in range(len_info):
            # Sequence
            item = QtWidgets.QTableWidgetItem(str(info[row][0]))
            self.table_m.setItem(row, 0, item)
            # Study acronym
            item = QtWidgets.QTableWidgetItem(info[row][1])
            self.table_m.setItem(row, 1, item)
            # Date and start time
            item = QtWidgets.QTableWidgetItem(info[row][2])
            self.table_m.setItem(row, 2, item)
            # Description
            item = QtWidgets.QTableWidgetItem(info[row][3])
            self.table_m.setItem(row, 3, item)
            # Filename
            item = QtWidgets.QTableWidgetItem(info[row][4])
            self.table_m.setItem(row, 4, item)
            self.file_list.append(info[row][4])
            # Xcen and Ycen
            fstring = '{:0.1f}'.format(info[row][5])
            item = QtWidgets.QTableWidgetItem(fstring)
            self.table_m.setItem(row, 5, item)
            fstring = '{:0.1f}'.format(info[row][6])
            item = QtWidgets.QTableWidgetItem(fstring)
            self.table_m.setItem(row, 6, item)

        # Any cells highlighted?
        self.table_m.cellClicked.connect(self.get_details)

    def get_details(self, row, column):
        """Provide details on selected cell."""
        self.info_detail.clear()
        info = self.fill_info(str(self.table_m.item(row, 4).text()))
        for line in info:
            self.info_detail.append(line)
        self.info_detail.verticalScrollBar().\
            setValue(self.info_detail.verticalScrollBar().minimum())

    def fill_info(self, file):
        """Retrieve useful info for a selected file."""
        info = []
        if len(self.d.eis_str) != 0:
            row, = [x for x in self.d.eis_str if x.filename == str(file)]
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
            temp = "".join(format(x, "<6.1f") for x in row.wave)
            info.append("{0:<20} {1}".format('wave', temp))
            temp = "".join(format(x, "<6.1f") for x in row.wavemin)
            info.append("{0:<20} {1}".format('wavemin', temp))
            temp = "".join(format(x, "<6.1f") for x in row.wavemax)
            info.append("{0:<20} {1}".format('wavemax', temp))
            info.append("{0:<20} {1}".format('width', row.width))
            info.append("{0:<20} {1}".format('obstitle', row.obstitle))
            info.append("{0:<20} {1}".format('obs_dec', row.obs_dec))
            info.append("{0:<20} {1}".format('sci_obj', row.sci_obj))
            info.append("{0:<20} {1}".format('slit_index', row.slit_index))
            info.append("{0:<20} {1}".format('scan_fm_nsteps', row.scan_fm_nsteps))
            info.append("{0:<20} {1}".format('scan_fm_stepsize',
                                             row.scan_fm_stepsize))
            info.append("{0:<20} {1}".format('nexp', row.nexp))
            info.append("{0:<20} {1}".format('exptime', row.exptime))
        return info

    def details(self):
        """Display detailed cat info."""
        title = QtWidgets.QLabel(self)
        title.setText('Details on Selected Item')
        self.grid.addWidget(title, self.gui_row, 0, 1, 6)
        self.gui_row += 1

        self.info_detail = QtWidgets.QTextEdit()
        if have_Qt4:
            font = QtWidgets.QFont("Courier New")
        else:
            font = QtGui.QFont("Courier New")
        self.info_detail.setFont(font)
        self.grid.addWidget(self.info_detail, self.gui_row, 0, 1, 6)
        #self.grid.setRowStretch(self.gui_row, 1)
        self.gui_row += 1
        self.info_detail.setReadOnly(True)

    def bottom(self):
        """Display bottom info."""
        download_list = QtWidgets.QPushButton('Download Files', self)
        download_list.setFixedWidth(self.default_button_width)
        self.grid.addWidget(download_list, self.gui_row, 0)
        download_list.clicked.connect(self.event_download_file_list)

        self.topdir_box = QtWidgets.QLineEdit(self)
        self.topdir_box.setFixedWidth(self.default_button_width)
        self.topdir_box.setText(self.default_topdir)
        self.grid.addWidget(self.topdir_box, self.gui_row, 1)

        self.radio = QtWidgets.QRadioButton("Use Date Tree")
        self.radio.setFixedWidth(self.default_button_width)
        self.grid.addWidget(self.radio, self.gui_row, 2)

        self.gui_row += 1

        self.save_list = QtWidgets.QPushButton('Save File List', self)
        self.save_list.setFixedWidth(self.default_button_width)
        self.grid.addWidget(self.save_list, self.gui_row, 0)
        self.save_list.clicked.connect(self.event_save_file_list)

        self.filename_box = QtWidgets.QLineEdit(self)
        self.filename_box.setFixedWidth(self.default_button_width)
        self.filename_box.setText(self.default_filename)
        self.grid.addWidget(self.filename_box, self.gui_row, 1)

    def event_download_file_list(self):
        if self.file_list is not None:
            datetree = self.radio.isChecked()
            topdir = self.topdir_box.text()

#            sys.stdout = OutLog(self.info_detail, sys.stdout)
#            sys.stderr = OutLog(self.info_detail, sys.stderr)

            o = download_hdf5_data(filename=self.file_list, datetree=datetree, local_top=topdir)

    def event_save_file_list(self):
        """Save the displayed list of files, one per line."""
        if self.file_list is not None:
            filename = self.filename_box.text()
            if filename != '':
                with open(filename, 'w') as fp:
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
        # option 2a - use file from SSW (configured environmentS)
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

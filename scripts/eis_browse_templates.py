#!/usr/bin/env python
__all__ = ['eis_browse_templates']

import sys
import os
import numpy as np
import shutil
from PyQt5 import Qt, QtCore, QtGui, QtWidgets
from PyQt5.QtGui import QIcon, QPixmap, QImage
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from eispac.templates import EISTemplateLocator

class Top(QWidget):

    def __init__(self, eis):
        super().__init__()
        self.eis = eis
        self.save_dir = os.path.join(os.getcwd(), 'eis_fit_templates') # 'eis_template_dir'
        self.font_default = QtGui.QFont("Courier", 12)
        self.font_small = QtGui.QFont("Courier", 10)
        self.imgNX = 1458
        self.imgNY = 614
        self.colNX = 200 # min width of the three columns at top the the gui

        self.initUI()

    def initUI(self):
        self.buttonQuit = QPushButton('Quit')
        self.buttonQuit.setFont(self.font_default)
        self.buttonQuit.clicked.connect(self.on_click_button_quit)
        self.buttonQuit.resize(self.colNX, self.frameGeometry().height())

        self.buttonSelect = QPushButton('Select header')
        self.buttonSelect.setFont(self.font_default)
        self.buttonSelect.clicked.connect(self.on_click_button_select)
        self.buttonSelect.resize(self.colNX, self.frameGeometry().height())

        self.buttonDir = QPushButton('Change save dir')
        self.buttonDir.setFont(self.font_default)
        self.buttonDir.clicked.connect(self.on_click_button_dir)
        self.buttonDir.resize(self.colNX, self.frameGeometry().height())

        self.buttonCopy = QPushButton('Copy template')
        self.buttonCopy.setFont(self.font_default)
        self.buttonCopy.clicked.connect(self.on_click_button_copy)
        self.buttonCopy.resize(self.colNX, self.buttonCopy.frameGeometry().height())

        self.fileLabel = QLabel()
        self.fileLabel.setFont(self.font_default)
        self.fileLabel.setAlignment(Qt.AlignCenter)
        self.fileLabel.resize(self.colNX, self.fileLabel.frameGeometry().height())
        self.setup_file_label()

        self.listTemplates = QListWidget()
        self.listTemplates.setFont(self.font_default)
        self.listTemplates.clicked.connect(self.on_click_list_templates)
        self.setup_list()

        self.textWindow = QTextEdit()
        self.textWindow.setFont(self.font_small)
        info = ('* Select an EIS HDF header file to see a list of relevent'
               +' templates\n\n'
               +'* Select a row in the table to display an example plot of'
               +' the template applied to representative solar spectra.\n\n'
               +'* Use the copy button to copy the selected template file to'
               +' the directory given below.\n\n'
               +'* Some templates contain multiple component Gaussians. Each'
               +' component is listed seperately, however you only need to copy'
               +' one file.')
        self.textWindow.append(info+'\n')
        self.textWindow.setReadOnly(True)
        self.textWindow.resize(self.colNX, self.textWindow.frameGeometry().height())

        self.savedirLabel = QLabel()
        self.savedirLabel.setFont(self.font_default)
        self.savedirLabel.setAlignment(Qt.AlignLeft)
        self.savedirLabel.resize(self.colNX, self.savedirLabel.frameGeometry().height())
        self.savedirLabel.setText('Save directory')

        self.savedirBox = QLineEdit(self)
        self.savedirBox.resize(self.colNX, self.savedirBox.frameGeometry().height())
        self.savedirBox.setText(self.save_dir)

        self.window = QLabel()
        buff = np.zeros((self.imgNX, self.imgNX, 3), dtype=np.int32)
        image = QImage(buff, self.imgNX, self.imgNY, QImage.Format_ARGB32)
        self.window.setPixmap(QPixmap(image))

        # TO-DO: Replace this mess of nested boxes with a proper widget grid
        left_box = QVBoxLayout()
        left_box.addWidget(self.buttonQuit)
        left_box.addWidget(self.buttonSelect)
        left_box.addWidget(self.buttonDir)
        left_box.addWidget(self.buttonCopy)
        left_box.addStretch()

        mid_box = QVBoxLayout()
        mid_box.addWidget(self.fileLabel)
        mid_box.addWidget(self.listTemplates)

        right_box = QVBoxLayout()
        right_box.addWidget(self.textWindow)
        right_box.addWidget(self.savedirLabel)
        right_box.addWidget(self.savedirBox)

        top_box = QHBoxLayout()
        top_box.addLayout(left_box)
        top_box.addLayout(mid_box)
        top_box.addLayout(right_box)

        bot_box = QVBoxLayout()
        bot_box.addWidget(self.window)

        vbox = QVBoxLayout()
        vbox.addLayout(top_box)
        vbox.addLayout(bot_box)

        # --- display the widget
        self.setLayout(vbox)
        self.setGeometry(50, 100, 1600, 900)
        self.setWindowTitle('Browse EIS Fitting Template')
        self.show()

    def setup_list(self):
        if self.eis.text_list is not None:
            self.listTemplates.clear()
            self.listTemplates.addItems(self.eis.text_list)
        self.listTemplates.setFont(self.font_default)
        count = self.listTemplates.count()
        if count == 0: count = 2
        if count > 15: count = 15
        nrows = self.listTemplates.sizeHintForRow(0)*count + \
            2*self.listTemplates.frameWidth() + 5
        # self.listTemplates.setFixedSize(400, nrows)
        self.listTemplates.resize(self.colNX, nrows)

    def setup_file_label(self):
        if self.eis.filename_head is not None:
            f = os.path.basename(self.eis.filename_head)
            f = f.split('.')[0]
            s = ('Selected header: '+f+'\n'
                +'template file,       window,   wmin,   wmax')
            self.fileLabel.setText(s)
        else:
            self.fileLabel.setText('No file has been selected')

    def on_click_button_quit(self):
        # --- quit the app
        QApplication.instance().quit()

    def on_click_button_select(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        eis_filename, _ = QFileDialog.getOpenFileName(self, 'Select a file',
                                                      filter='*.head.h5',
                                                      options=options)
        if os.path.isfile(eis_filename):
            self.eis = EISTemplateLocator(eis_filename, ignore_local=True)
            self.setup_list()
            self.setup_file_label()
            self.set_blank_image()

    def on_click_button_dir(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        options |= QFileDialog.ShowDirsOnly
        new_savedir = QFileDialog.getExistingDirectory(self, 'Select a directory',
                                                       options=options)
        if os.path.isdir(new_savedir):
            self.save_dir = new_savedir
            self.savedirBox.setText(self.save_dir)

    def on_click_button_copy(self):
        item = self.listTemplates.currentItem()
        index = self.listTemplates.indexFromItem(item)
        n = index.row()
        template_file = self.eis.template_list[n][0]

        self.save_dir = self.savedirBox.text()
        if not os.path.isdir(self.save_dir):
            os.mkdir(self.save_dir)
            s = f'Created {self.save_dir}'
            print(s)
            self.textWindow.append(s)

        shutil.copy2(template_file, self.save_dir)

        s = f'Copied {os.path.basename(template_file)} to {self.save_dir}'
        print(s)
        self.textWindow.append(s)

    def on_click_list_templates(self):
        item = self.listTemplates.currentItem()
        index = self.listTemplates.indexFromItem(item)
        n = index.row()
        image_file = self.eis.template_list[n][0].replace('.h5', '.jpg')
        self.display_image(image_file)

    def display_image(self, image_file):
        if os.path.isfile(image_file):
            pixmap = QPixmap(image_file)
            pixmap = pixmap.scaled(self.imgNX, self.imgNY)
            self.window.setPixmap(pixmap)
        else:
            s = f'Error: cannot find or display {image_file}'
            self.textWindow.append(s)

    def set_blank_image(self):
        buff = np.zeros((self.imgNX, self.imgNX, 3), dtype=np.int32)
        image = QImage(buff, self.imgNX, self.imgNY, QImage.Format_ARGB32)
        self.window.setPixmap(QPixmap(image))

def eis_browse_templates():
    # check the input
    if len(sys.argv) > 1:
        eis_filename = sys.argv[1]
    else:
        eis_filename = None

    # create the object that actually finds the templates
    eis = EISTemplateLocator(eis_filename, ignore_local=True)

    app = QApplication(sys.argv)
    ex = Top(eis)
    sys.exit(app.exec_())

if __name__ == '__main__':
    eis_browse_templates()

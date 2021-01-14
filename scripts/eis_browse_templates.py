#!/Users/hpw/Desktop/projects/eispac/software/local-anaconda3/bin/python
__all__ = ['eis_browse_templates']

import sys
import os
import numpy as np
import shutil
from PyQt5 import Qt, QtCore, QtGui, QtWidgets
from PyQt5.QtGui import QIcon, QPixmap, QImage
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from eispac.templates import eis_find_templates

class Top(QWidget):

    def __init__(self, eis):
        super().__init__()
        self.eis = eis
        self.default_dir = 'eis_template_dir'
        self.font_default = QtGui.QFont("Courier", 12)
        self.font_small = QtGui.QFont("Courier", 10)
        self.winNX = 1458
        self.winNY = 614
        self.leftNX = 450

        self.initUI()

    def initUI(self):
        self.buttonQuit = QPushButton('Quit')
        self.buttonQuit.setFont(self.font_default)
        self.buttonQuit.clicked.connect(self.on_click_button_quit)
        self.buttonQuit.resize(self.leftNX, self.frameGeometry().height())

        self.buttonSelect = QPushButton('Select a header file')
        self.buttonSelect.setFont(self.font_default)
        self.buttonSelect.clicked.connect(self.on_click_button_select)
        self.buttonSelect.resize(self.leftNX, self.frameGeometry().height())

        self.fileLabel = QLabel()
        self.fileLabel.setFont(self.font_default)
        self.fileLabel.setAlignment(Qt.AlignCenter)
        self.fileLabel.resize(self.leftNX, self.fileLabel.frameGeometry().height())
        self.setup_file_label()

        self.listTemplates = QListWidget()
        self.listTemplates.setFont(self.font_default)
        self.listTemplates.clicked.connect(self.on_click_list_templates)
        self.setup_list()

        self.buttonCopy = QPushButton('Copy Template to Local Dir')
        self.buttonCopy.setFont(self.font_default)
        self.buttonCopy.clicked.connect(self.on_click_button_copy)
        self.buttonCopy.resize(self.leftNX, self.buttonCopy.frameGeometry().height())

        self.textWindow = QTextEdit()
        self.textWindow.setFont(self.font_small)
        info = ('* Select an EIS HDF header file.\n\n'
            +'* Templates relevant to that file will be listed.\n\n'
            +'* Make a selection to display the template applied to some'
            +' represenative solar spectra.\n\n'
            +'* Use the copy button to copy the template file to a local'
            +' directory.\n\n'
            +'* You only need to copy one file for the template of a'
            ' multi-component fit. All components are listed separately.')
        self.textWindow.append(info+'\n')
        self.textWindow.setReadOnly(True)
        self.textWindow.resize(self.leftNX, self.textWindow.frameGeometry().height())

        self.window = QLabel()
        buff = np.zeros((self.winNX, self.winNX, 3), dtype=np.int16)
        image = QImage(buff, self.winNX, self.winNY, QImage.Format_ARGB32)
        self.window.setPixmap(QPixmap(image))

        vbox1 = QVBoxLayout()
        vbox1.addWidget(self.buttonQuit)
        vbox1.addWidget(self.buttonSelect)
        vbox1.addWidget(self.fileLabel)
        vbox1.addWidget(self.listTemplates)
        vbox1.addWidget(self.buttonCopy)
        vbox1.addWidget(self.textWindow)
        vbox1.addStretch()

        vbox2 = QVBoxLayout()
        vbox2.addWidget(self.window)

        hbox = QHBoxLayout()
        hbox.addLayout(vbox1)
        hbox.addLayout(vbox2)

        self.setLayout(hbox)

        # --- display the widget
        self.setWindowTitle('Select EIS Fitting Template')
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
        self.listTemplates.resize(self.leftNX, nrows)


    def setup_file_label(self):
        if self.eis.filename_head is not None:
            f = os.path.basename(self.eis.filename_head)
            f = f.split('.')[0]
            s = ('Selected header: '+f+'\n'
                +'--- Available Templates ---\n'
                +'filename,   window num,   wmin,   wmax')
            self.fileLabel.setText(s)
        else:
            self.fileLabel.setText('No file has been selected')

    def on_click_button_quit(self):
        # --- quit the app
        QApplication.instance().quit()

    def on_click_button_copy(self):
        item = self.listTemplates.currentItem()
        index = self.listTemplates.indexFromItem(item)
        n = index.row()
        template_file = self.eis.template_list[n][0]

        if not os.path.isdir(self.default_dir):
            os.mkdir(self.default_dir)
            s = f'created {self.default_dir}'
            print(s)
            self.textWindow.append(s)

        shutil.copy2(template_file, self.default_dir)

        s = f'copied {os.path.basename(template_file)} to {self.default_dir}'
        print(s)
        self.textWindow.append(s)

    def on_click_button_select(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        eis_filename, _ = QFileDialog.getOpenFileName(self, 'Select a file',
                                                      filter='*.head.h5',
                                                      options=options)
        if os.path.isfile(eis_filename):
            self.eis = eis_find_templates(eis_filename, ignore_local=True)
            self.setup_list()
            self.setup_file_label()
            self.set_blank_image()

    def on_click_list_templates(self):
        item = self.listTemplates.currentItem()
        index = self.listTemplates.indexFromItem(item)
        n = index.row()
        image_file = self.eis.template_list[n][0].replace('.h5', '.jpg')
        self.display_image(image_file)

    def display_image(self, image_file):
        if os.path.isfile(image_file):
            pixmap = QPixmap(image_file)
            pixmap = pixmap.scaled(self.winNX, self.winNY)
            self.window.setPixmap(pixmap)

    def set_blank_image(self):
        buff = np.zeros((self.winNX, self.winNX, 3), dtype=np.int16)
        image = QImage(buff, self.winNX, self.winNY, QImage.Format_ARGB32)
        self.window.setPixmap(QPixmap(image))

def eis_browse_templates():
    # check the input
    if len(sys.argv) > 1:
        eis_filename = sys.argv[1]
    else:
        eis_filename = None

    # create the object that actually finds the templates
    eis = eis_find_templates(eis_filename, ignore_local=True)

    app = QApplication(sys.argv)
    ex = Top(eis)
    sys.exit(app.exec_())

if __name__ == '__main__':
    eis_browse_templates()

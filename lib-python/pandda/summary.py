import os, sys

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtWebKit import *

class PanddaHtmlWidget(QWidget):
    def __init__(self, name, content):
        super(PanddaHtmlWidget, self).__init__()
        # Widget Name
        self.name = name
        self.content = content

        # Add a VBox to the panddaWidget
        self.layout = QVBoxLayout(self)
        self.view = None

#        # Add refresh button
#        self.refresh_button = QPushButton('Refresh', self)
#        self.refresh_button.clicked.connect(self.refresh)
#        # Add the button to the vbox
#        self.layout.addWidget(self.refresh_button)

        # Add the html
        self.build_html()

    def build_html(self):
        # Delete old view
        if self.view: self.view.deleteLater()
        # Add a view to the panddaWidget
        self.view = QWebView(self)
        self.layout.addWidget(self.view)
        # Add html to the view
        if os.path.exists(self.content):
            with open(self.content, 'r') as fh:
                text = fh.read()
            dir = os.path.dirname(os.path.abspath(self.content))
            # Dynamically change the html paths to full paths to allow rendering
            text = text.replace('src="../', 'src="file://'+dir+'/../')
            text = text.replace('src="./', 'src="file://'+dir+'/')
            # Set the text in the html window
            self.view.setHtml("<html><head></head><body></body></html>")
            self.view.setHtml(text)
        else:
            self.view.setHtml('File does not yet exist: {}'.format(self.content))
        self.view.show()

    def refresh(self):
        print 'Refreshing html'
        self.build_html()

class PanddaTabBar(QTabWidget):
    def __init__(self, name=None, tabs=None):
        super(PanddaTabBar, self).__init__()
        self.name = name
        # Create Tabs
        self.tabs = []
        if tabs is not None:
            for t in tabs:
                self.add_tab(t)
    def add_tab(self, widget):
        self.tabs.append(widget)
        self.addTab(widget, widget.name)

class PanddaWindow(PanddaTabBar):
    def __init__(self, tabs=None, title='PanDDA HTML Summaries'):
        super(PanddaWindow, self).__init__(name=None, tabs=tabs)
        self.setWindowTitle(title)


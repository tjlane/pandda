import os, sys

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtWebKit import *

class PanddaHtmlWidget(QWidget):
    def __init__(self, name, content):
        super(PanddaHtmlWidget, self).__init__()

        self.name = name
        self.content = content
        self.layout = QVBoxLayout(self)
        self.view = None
        self._built = False

    def is_built(self):
        return self._built

    def build_html(self):
        print '> Loading HTML: {} ({})'.format(self.name, self.content)
        if self.view: self.view.deleteLater()
        self.view = QWebView(self)
        self.layout.addWidget(self.view)
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
        self._built = True
        self.view.show()

    def refresh(self):
        print 'Refreshing html'
        self.build_html()

class PanddaTabBar(QTabWidget):
    def __init__(self, name=None, tabs=None):
        super(PanddaTabBar, self).__init__()

        self.name = name
        self.tabs = []

        if tabs is not None:
            for t in tabs:
                self.add_tab(t)

        self.currentChanged.connect(self.load_tab)

    def add_tab(self, widget):
        self.tabs.append(widget)
        self.addTab(widget, widget.name)

    def load_tab(self, i):
        t = self.tabs[i]
        if not t.is_built():
            t.build_html()

class PanddaWindow(QTabWidget):
    def __init__(self, tabs=None, title='PanDDA HTML Summaries'):
        super(PanddaWindow, self).__init__()

        self.setWindowTitle(title)

        self.tabs = []

        if tabs is not None:
            for t in tabs:
                self.add_tab(t)

    def add_tab(self, widget):
        if isinstance(widget, PanddaHtmlWidget):
            widget.build_html()
        self.tabs.append(widget)
        self.addTab(widget, widget.name)


import os, sys

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtWebKit import *

class panddaHtmlWidget(QWidget):
    def __init__(self, name, content):
        super(panddaHtmlWidget, self).__init__()
        # Widget Name
        self.name = name
        # Add a view to the panddaWidget
        view = QWebView(self)
        # Add a VBox to the panddaWidget
        layout = QVBoxLayout(self)
        # Add the view to the vbox
        layout.addWidget(view)
        # Add html to the view
        if os.path.exists(content):
            with open(content, 'r') as fh:
                text = fh.read()
            dir = os.path.dirname(os.path.abspath(content))
            # Dynamically change the html paths to full paths to allow rendering
            text = text.replace('src="../', 'src="file://'+dir+'/../')
            text = text.replace('src="./', 'src="file://'+dir+'/')
            # Set the text in the html window
            view.setHtml(text)
        else:
            view.setHtml('File does not yet exist: {}'.format(content))

class panddaWindow(QTabWidget):
    def __init__(self, tabs=None):
        super(panddaWindow, self).__init__()
        self.setWindowTitle('PANDDA GUI')
        # Create Tabs
        self.tabs = []
        if tabs is not None:
            for t in tabs:
                self.add_tab(t)
    def add_tab(self, widget):
        self.tabs.append(widget)
        self.addTab(widget, widget.name)

#class Render(QWebPage):
#  def __init__(self, url):
#    self.app = QApplication(sys.argv)
#    QWebPage.__init__(self)
#    self.loadFinished.connect(self._loadFinished)
#    self.mainFrame().load(QUrl(url))
#    self.app.exec_()
#
#  def _loadFinished(self, result):
#    self.frame = self.mainFrame()
#    self.app.quit()

def test_1(test_dir):
    app = QApplication(sys.argv)
    window = panddaWindow()
    window.add_tab(panddaHtmlWidget(name='RESULTS', content=os.path.join(test_dir,'./pandda_results_example_basic.html')))
    window.show()
    app.exec_()

#def test_2(test_dir):
#    url = path2url(os.path.join(test_dir,'pandda_results_example.html'))
#    print url
#    html = Render(url).frame.toHtml()
#    print html
#
#    app = QApplication(sys.argv)
#    window = panddaWindow()
#    widget = panddaHtmlWidget(name='RESULTS', content=html)
#    window.add_tab(widget)
#    window.show()
#    app.exec_()

if __name__ == "__main__":
    test_dir = os.path.join(os.path.dirname(__file__), 'test_files')
    test_1(test_dir)
    test_2(test_dir)

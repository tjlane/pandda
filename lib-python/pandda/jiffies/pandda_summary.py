import os, sys
import threading

from PyQt4.QtGui import QApplication
from pandda.summary import panddaWindow, panddaHtmlWidget

#######################################

blank_arg_prepend = None

master_phil = None

HTML_1 = './analyses/html_summaries/pandda_initial.html'
HTML_2 = './analyses/html_summaries/pandda_analyse.html'
HTML_3 = './analyses/html_summaries/pandda_inspect.html'

#######################################

def show_summary():
    app = QApplication(sys.argv)
    window = panddaWindow()
    window.add_tab(panddaHtmlWidget(name='Dataset Summary', content=HTML_1))
    window.add_tab(panddaHtmlWidget(name='Results Summary', content=HTML_2))
    window.add_tab(panddaHtmlWidget(name='Inspect Summary', content=HTML_3))
    # Choose which window is initally open
    if os.path.exists(HTML_3):
        window.setCurrentIndex(2)
    elif os.path.exists(HTML_2):
        window.setCurrentIndex(1)
    window.show()
    app.exec_()

#######################################

def run():
    t = threading.Thread(target=show_summary, args=())
    t.daemon = True
    t.start()
    t.join()

#######################################

if __name__=='__main__':
    run()

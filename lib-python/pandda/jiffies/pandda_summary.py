import os, sys

from PyQt4.QtGui import QApplication
from pandda.summary import panddaWindow, panddaHtmlWidget

#######################################

try:
    sys.path.append(os.path.join(os.environ['PANDDA_TOP'], 'lib-python'))
except:
    pass

#######################################

blank_arg_prepend = None

master_phil = None

#######################################

def show_summary():
    app = QApplication(sys.argv)
    window = panddaWindow()
    window.add_tab(panddaHtmlWidget(name='Dataset Summary', content='./results_summaries/pandda_initial.html'))
    window.add_tab(panddaHtmlWidget(name='Results Summary', content='./results_summaries/pandda_analyse.html'))
    window.add_tab(panddaHtmlWidget(name='Inspect Summary', content='./results_summaries/pandda_inspect.html'))
    # Choose which window is initally open
    if   os.path.exists('./results_summaries/pandda_inspect.html'):
        window.setCurrentIndex(2)
    elif os.path.exists('./results_summaries/pandda_analyse.html'):
        window.setCurrentIndex(1)
    window.show()
    app.exec_()

#######################################

def run(run=None, master_phil=None, args=None, blank_arg_prepend=None):
    show_summary()

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)

import os, sys, glob
import threading

from PyQt4.QtGui import QApplication
from pandda.summary import PanddaWindow, PanddaHtmlWidget, PanddaTabBar

#######################################

blank_arg_prepend = None

master_phil = None

INIT_HTML = './analyses/html_summaries/pandda_initial.html'
MAP_HTMLS = sorted(glob.glob('./analyses/html_summaries/pandda_map_*.html'))
ANAL_HTML = './analyses/html_summaries/pandda_analyse.html'
INSP_HTML = './analyses/html_summaries/pandda_inspect.html'

#######################################

def show_summary():
    app = QApplication(sys.argv)
    window = PanddaWindow(title='PanDDA HTML Summaries: {}'.format(os.getcwd()))
    # ----------------------------------
    window.add_tab(PanddaHtmlWidget(name='Dataset Summary', content=INIT_HTML))
    # ----------------------------------
    multi_tab = PanddaTabBar(name='Map Analysis Summaries')
    for f in MAP_HTMLS:
        name = f[f.find('pandda_map_'):].replace('pandda_map_','').replace('.html','')
        multi_tab.add_tab(PanddaHtmlWidget(name=name, content=f))
    window.add_tab(multi_tab)
    # ----------------------------------
    window.add_tab(PanddaHtmlWidget(name='Results Summary', content=ANAL_HTML))
    # ----------------------------------
    window.add_tab(PanddaHtmlWidget(name='Inspect Summary', content=INSP_HTML))
    # ----------------------------------
    # Choose which window is initally open
    if os.path.exists(INSP_HTML):
        window.setCurrentIndex(3)
    elif os.path.exists(ANAL_HTML):
        window.setCurrentIndex(2)
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

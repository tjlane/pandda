import os, sys, glob
import threading

from PyQt4.QtGui import QApplication
from pandda.summary import PanddaWindow, PanddaHtmlWidget, PanddaTabBar

#######################################

blank_arg_prepend = None

master_phil = None

INIT_HTML = './analyses/html_summaries/pandda_initial.html'
MAP_HTMLS = sorted(glob.glob('./analyses/html_summaries/pandda_map_*.html'))
DST_HTMLS = sorted(glob.glob('./processed_datasets/*/html/*.html'))
ANAL_HTML = './analyses/html_summaries/pandda_analyse.html'
INSP_HTML = './analyses/html_summaries/pandda_inspect.html'

#######################################

def show_summary():
    app = QApplication(sys.argv)
    window = PanddaWindow(title='PanDDA HTML Summaries: {}'.format(os.getcwd()))
    # ----------------------------------
    print '> Creating Tab: '+INIT_HTML
    INIT_WIDG = PanddaHtmlWidget(name='Initial Summary', content=INIT_HTML)
    window.add_tab(INIT_WIDG)
    # ----------------------------------
    if DST_HTMLS:
        multi_tab = PanddaTabBar(name='Individual Dataset Summaries')
        for f in DST_HTMLS:
            print '> Creating Tab: '+f
            name = os.path.splitext(os.path.basename(f))[0]
            multi_tab.add_tab(PanddaHtmlWidget(name=name, content=f))
        window.add_tab(multi_tab)
    # ----------------------------------
    if MAP_HTMLS:
        multi_tab = PanddaTabBar(name='Map Analysis Summaries')
        for f in MAP_HTMLS:
            print '> Creating Tab: '+f
            name = os.path.splitext(os.path.basename(f))[0].replace('pandda_map_','')
            multi_tab.add_tab(PanddaHtmlWidget(name=name, content=f))
        window.add_tab(multi_tab)
    # ----------------------------------
    print '> Creating Tab: '+ANAL_HTML
    ANAL_WIDG = PanddaHtmlWidget(name='Results Summary', content=ANAL_HTML)
    window.add_tab(ANAL_WIDG)
    # ----------------------------------
    print '> Creating Tab: '+INSP_HTML
    INSP_WIDG = PanddaHtmlWidget(name='Inspect Summary', content=INSP_HTML)
    window.add_tab(INSP_WIDG)

    # ----------------------------------
    # Choose which window is initally open
    if os.path.exists(INSP_HTML):
        window.setCurrentWidget(INSP_WIDG)
    elif os.path.exists(ANAL_HTML):
        window.setCurrentWidget(ANAL_WIDG)
    else:
        window.setCurrentWidget(INIT_WIDG)
    window.show()
    app.exec_()

#######################################

def run():
    show_summary()

#######################################

if __name__=='__main__':
    run()

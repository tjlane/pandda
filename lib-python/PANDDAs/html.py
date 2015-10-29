
import os, sys
import urlparse, urllib
import jinja2

PANDDA_HTML_ENV = jinja2.Environment(loader=jinja2.PackageLoader('PANDDAs', 'templates'))

def path2url(path):
    return urlparse.urljoin('file:', urllib.pathname2url(path))


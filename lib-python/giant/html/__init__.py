from future import standard_library
standard_library.install_aliases()

import giant.logs as lg
logger = lg.getLogger(__name__)

import os
import jinja2
import urllib.parse, urllib.request, urllib.parse, urllib.error

from . import divs
from . import summary


HTML_ENV = jinja2.Environment(
    loader = jinja2.PackageLoader('giant', 'templates')
    )

def path2url(path):
    return urllib.parse.urljoin('file:', urllib.request.pathname2url(path))

def png2base64str(path):
    """Convert a png for embedding into an html page"""
    contents = open(path, 'rb').read().encode('base64').replace('\n', '')
    return contents

def png2base64src(path):
    return 'data:image/png;base64,{}'.format(png2base64str(path))

def png2base64src_maybe(path, print_on_missing=False):
    if (path is not None) and os.path.exists(path):
        return png2base64src(path)
    else:
        if print_on_missing:
            logger('missing file: {}'.format(path))
        return 'none'

class AutoDivNamer(object):
    def __init__(self, start_index=1):
        self.index=start_index

    def format(self, index):
        return 'div{}'.format(index)

    def current(self):
        return self.format(self.index)

    def __next__(self):
        self.index += 1
        return self.current()

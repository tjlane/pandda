from future import standard_library
standard_library.install_aliases()

import giant.logs as lg
logger = lg.getLogger(__name__)

import os
import base64
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
    with open(path, 'rb') as fh:
        contents = fh.read()
    b64_str = base64.b64encode(contents)
    return b64_str

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

    def next(self):
        self.index += 1
        return self.current()

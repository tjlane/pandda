
import jinja2
import urlparse, urllib

BAMBOO_HTML_ENV = jinja2.Environment(loader=jinja2.PackageLoader('bamboo', 'templates'))

def path2url(path):
    return urlparse.urljoin('file:', urllib.pathname2url(path))


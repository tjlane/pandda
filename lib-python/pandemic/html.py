
import jinja2

class _AutoDivNamer(object):
    def __init__(self, start_index=1):
        self.index=start_index

    def format(self, index):
        return 'div{}'.format(index)

    def current(self):
        return self.format(self.index)

    def next(self):
        self.index += 1
        return self.current()

PANDEMIC_HTML_ENV = jinja2.Environment(loader=jinja2.PackageLoader('pandemic', 'templates'))

PANDEMIC_HTML_ENV.globals['divnamer'] = _AutoDivNamer

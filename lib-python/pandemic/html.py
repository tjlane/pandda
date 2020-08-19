
import jinja2

from giant.html import AutoDivNamer

class PandemicAutoDivNamer(AutoDivNamer):
    pass

HTML_ENV = jinja2.Environment(
    loader=jinja2.ChoiceLoader([
        jinja2.PackageLoader('pandemic', 'templates'),
        jinja2.PackageLoader('giant', 'templates'),
        ]),
    extensions=['jinja2.ext.do'],
    trim_blocks=True,
    lstrip_blocks=True,
    )

HTML_ENV.globals['divnamer'] = PandemicAutoDivNamer

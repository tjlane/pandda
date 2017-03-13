import os, sys, datetime


class Bar(object):


    def __init__(self, width=40, body='-', head='>>>'):
        self.width = width
        self.body  = body
        self.head  = head

    def __call__(self):
        width = self.width - len(self.head)
        if width < 1: width = 0
        return self.body*width + self.head


class Heading(object):


    def __init__(self, width=80, spacer='#', decorator=' <~~~> ', edge_width=3):
        """Spacer for printing neatly formatted headings"""

        assert width > 10
        self._decorator     = decorator
        self._spacer        = spacer
        self._width         = width
        self._content_width = width-(2*edge_width)
        self._side_width    = edge_width

        self._side_padding  = self._side_width*spacer

        self._blank_line = ''
        self._text_line  = self._side_padding + '{{:^{}}}'.format(self._content_width) + self._side_padding
        self._inner_line = self._text_line.format('')
        self._outer_line = self._decorator.center(self._width, self._spacer)

    def __call__(self, text='', spacer=False, blank=False):
        out = []
        out.append(self._text_line.format(text))
        if spacer: out = [self._inner_line] + out + [self._inner_line]
        out = [self._outer_line] + out + [self._outer_line]
        if blank: out = [self._blank_line] + out + [self._blank_line]
        return '\n'.join(out)


class Log(object):


    def __init__(self, log_file=None, stdout=sys.stdout, verbose=False):
        """Log Object for writing to logs...?"""
        # File that the log will be written to
        self.log_file = log_file
        self.stdout = stdout
        self.verbose = verbose

        # Set default heading
        self.set_bar_and_heading(bar=Bar(), heading=Heading())

    def __call__(self, message, show=False, hide=False):
        """Log message to file, and mirror to stdout if verbose or force_print (hide overrules show)"""
        message = str(message)
        if (not hide) and (show or self.verbose):
            self.show(message)
        self.write(message=message+'\n', mode='a')

    def heading(self, message, spacer=False, blank=True):
        """Print message as a heading/divider"""
        text = self._heading(text=str(message), spacer=spacer, blank=blank)
        self.show(text)
        self.write(text)

    def bar(self, blank_before=False, blank_after=False):
        """Print divider bar"""
        text = '\n'*blank_before + self._bar() + '\n'*blank_after
        self.show(text)
        self.write(text)

    def set_bar_and_heading(self, bar=None, heading=None):
        if bar is not None:
            assert isinstance(bar, Bar)
            self._bar = bar
        if heading is not None:
            assert isinstance(heading, Heading)
            self._heading = heading
        return self

    def show(self, message):
        print(message)

    def write(self, message, mode='a'):
        if not self.log_file: return
        message = message.replace('\r','')
        with open(self.log_file, mode) as fh: fh.write(message)

    def read_all(self):
        return open(self.log_file, 'r').read()


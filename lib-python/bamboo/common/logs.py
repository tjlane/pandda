import os, sys, datetime

from bamboo.common import ListStream


class Bar(object):


    def __init__(self, width=40, body='-', head='>>>'):
        self.width = width
        self.body  = body
        self.head  = head

    def __call__(self, width=None, body=None, head=None):
        # Allow defaluts to be overridden
        if width is None: width = self.width
        if body is None:  body  = self.body
        if head is None:  head  = self.head
        # Modify width to allow for head
        width = width - len(head)
        if width < 1: width = 0
        return body*width + head


class Heading(object):


    def __init__(self, width=100, spacer='#', decorator=' <~~~> ', side_width=3):
        """Spacer for printing neatly formatted headings"""

        assert width > 10
        self._decorator     = decorator
        self._spacer        = spacer
        self._width         = width
        self._side_width    = side_width
        self._content_width = width-(2*self._side_width)

        self._side_padding  = self._side_width*spacer

    def __call__(self, text='', spacer=False, blank=False):

        # Allow for need to override defaults
        content_width = max(len(text)+2, self._content_width)
        total_width = content_width + (2*self._side_width)

        out = [self._text_line(text=text, text_width=content_width)]
        if spacer:
            out.insert(0, self._inner_line(text_width=content_width))
            out.append(out[0])
        out.insert(0, self._outer_line(width=total_width))
        out.append(out[0])
        if blank:
            out.insert(0, self._blank_line())
            out.append(out[0])
        return '\n'.join(out)

    def _blank_line(self):
        return ''

    def _text_line(self, text, text_width):
        return self._side_padding + text.center(text_width, ' ') + self._side_padding

    def _inner_line(self, text_width):
        return self._text_line(text='', text_width=text_width)

    def _outer_line(self, width):
        return self._decorator.center(width, self._spacer)


class ScreenLogger(object):


    def __init__(self, stdout=True):
        """Log Object for writing to logs...?"""
        # File that the log will be written to
        self.outputs = []

        # Set default heading
        self.set_bar_and_heading(bar        = Bar(),
                                 heading    = Heading(spacer='#', decorator=' <~~~> '),
                                 subheading = Heading(spacer='-', decorator=' *** '))

        if stdout is True:
            self.add_output(StdoutStream())

    def __call__(self, message):
        """Log message to file, and mirror to stdout if not silent and verbose or show is True"""
        message = str(message)+'\n'
        for o in self.outputs:
            o.write(message)

    def add_output(self, output):
        self.outputs.append(output)
        return self

    def heading(self, message, spacer=False, blank=True):
        """Print message as a heading/divider"""
        text = self._heading(text=str(message), spacer=spacer, blank=blank)
        self(text)

    def subheading(self, message, spacer=False, blank=True):
        """Print message as a heading/divider"""
        text = self._subheading(text=str(message), spacer=spacer, blank=blank)
        self(text)

    def bar(self, blank_before=False, blank_after=False, width=None, body=None, head=None):
        """Print divider bar"""
        text = '\n'*blank_before + self._bar(width=width, body=body, head=head) + '\n'*blank_after
        self(text)

    def set_bar_and_heading(self, bar=None, heading=None, subheading=None):
        if bar is not None:
            assert isinstance(bar, Bar)
            self._bar = bar
        if heading is not None:
            assert isinstance(heading, Heading)
            self._heading = heading
        if subheading is not None:
            assert isinstance(subheading, Heading)
            self._subheading = subheading
        return self

    def toggle(self, status=None):
        for o in self.outputs:
            if hasattr(o, 'toggle'):
                o.toggle(status)
        return self


class Log(ScreenLogger):


    def __init__(self, log_file=None, stdout=True):
        super(Log, self).__init__(stdout=stdout)
        if log_file is not None:
            self.add_output(LogFile(path=log_file))

    def log_file(self):
        """Returns the first LogFile object it contains (or None)"""
        for o in self.outputs:
            if isinstance(o, LogFile) or isinstance(o,LogStream):
                return o
        return None


class StdoutStream(object):


    def __init__(self, status=1):
        self._set_status(status)

    def _set_status(self, status):
        assert status in [0,1,False,True]
        self._status = int(status)

    def toggle(self, status=None):
        if status is None:
            self._set_status(1-self._status)
        else:
            self._set_status(status)

    def write(self, message):
        if self._status:
            sys.stdout.write(message)


class LogFile(object):


    def __init__(self, path):
        assert path is not None
        self.path = path

    def __str__(self):
        return self.path

    def as_log_stream(self):
        return LogStream(path=self.path)

    def write(self, message, mode='a'):
        message = message.replace('\r','')
        with open(self.path, mode) as fh:
            fh.write(message)

    def delete(self):
        """Delete the current log file"""
        if self.path and os.path.exists(self.path):
            os.remove(self.path)
        return self

    def read_all(self):
        if os.path.exists(self.path):
            return open(self.path, 'r').read()
        else:
            return None


class LogStream(ListStream):


    def __init__(self, path=None):
        super(LogStream, self).__init__()
        self.path = path

    def __str__(self):
        return self.path

    def as_log_file(self):
        return LogFile(path=self.path)

    def write_to_log(self, mode='a', clear_data=True):
        """Write out the stored contents (and delete contents if clear_data is True)"""
        assert self.path is not None
        self.as_log_file().write(self.format(), mode=mode)
        if clear_data is True:
            self.data = []
        return self



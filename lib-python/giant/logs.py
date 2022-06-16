import sys
import logging as lg

INFO = lg.INFO
WARNING = lg.WARNING
DEBUG = lg.DEBUG


class Bar(object):

    width = 40
    body = '-'
    head = '>>>'

    def __init__(self):
        self.body_width = max(0, self.width - len(self.head))
        self._bar = self.body*self.body_width + self.head

    def __call__(self, blank_before=False, blank_after=False):
        return '\n'*blank_before + self._bar + '\n'*blank_after


class Heading(object):

    width = 100
    spacer = '#'
    decorator = ' <~~~> '
    side_width = 3

    def __init__(self):
        self.content_width = max(0, self.width - (2*self.side_width))
        self.side_padding = self.side_width * self.spacer

    def __call__(self, text, spacer=False, blank=False):
        text = str(text)
        actual_width = max(len(text)+2, self.content_width)

        lines = [
            self.fmt(txt='',   width=actual_width, fill=self.spacer),
            self.fmt(txt=text, width=actual_width, fill=' '),
            self.fmt(txt='',   width=actual_width, fill=self.spacer),
        ]

        if spacer is True:
            s = self.fmt(txt='', width=actual_width, fill=' ')
            lines.insert(1, s)
            lines.insert(3, s)

        if blank is True:
            lines.insert(0, '\n')
            lines.append('\n')

        return '\n'.join(lines)

    def fmt(self, txt, width, fill=' '):
        return self.side_padding + txt.center(width, fill) + self.side_padding


class SubHeading(Heading):

    spacer = '-'
    decorator = ' ** '
    side_width = 3


class ListStream(object):
    """Replaces a file object for functions that print to screen - writes to list of lines instead"""

    _bar = Bar()
    _heading = Heading()
    _subheading = SubHeading()

    def __init__(self):
        self.data = []

    def __call__(self, s):
        self.write(s)

    def __str__(self):
        return self.format()

    # def __repr__(self):
    #     return self.format()

    def __iter__(self):
        return iter(self.data)

    def format(self):
        return ''.join(self.data)

    def write(self, s):
        self.data.append(s)

    def heading(self, message, spacer=False, blank=True):
        self.write(self._heading(text=message, spacer=spacer, blank=blank))

    def subheading(self, message, spacer=False, blank=True):
        self.write(self._subheading(text=message, spacer=spacer, blank=blank))

    def bar(self, blank_before=False, blank_after=False):
        self.write(self._bar(blank_before=blank_before, blank_after=blank_after))


class LoggerWithHeadings(lg.Logger):

    _bar = Bar()
    _heading = Heading()
    _subheading = SubHeading()

    def __call__(self, msg, *args, **kwargs):
        self.info(msg, *args, **kwargs)

    def heading(self, message, spacer=False, blank=True, level=lg.INFO):
        self.log(level=level, msg=self._heading(text=message, spacer=spacer, blank=blank))

    def subheading(self, message, spacer=False, blank=True, level=lg.INFO):
        self.log(level=level, msg=self._subheading(text=message, spacer=spacer, blank=blank))

    def bar(self, blank_before=False, blank_after=False, level=lg.INFO):
        self.log(level=level, msg=self._bar(blank_before=blank_before, blank_after=blank_after))


class ListHandler(lg.Handler): # Inherit from logging.Handler

    def __init__(self):
        # run the regular Handler __init__
        lg.Handler.__init__(self)
        # Custom attributes
        self.log_list = []

    def emit(self, record):
        # record.message is the log message
        self.log_list.append(record.msg)

    def size(self):
        return len(self.log_list)

    def list(self, start_i=0):
        return iter(self.log_list[start_i:])


class WarningListHandler(ListHandler):

    def __init__(self, name='warnings'):
        # run the regular ListHandler __init__
        ListHandler.__init__(self)
        # Custom attributes
        self.set_name(name)
        self.counter = 0

    def report_new(self, logger_name):
        assert logger_name is not None
        l = lg.getLogger(logger_name)
        # Get counter for next unreported message
        i = self.counter
        # Check if there are no new messages
        if i == self.size():
            return
        # Iterate through messages and report
        n = self.size() - i
        l.subheading('{} new warnings generated'.format(n))
        for w in self.list(start_i=i):
            l.bar(True, True)
            l(w)
        l.bar(True, True)
        # Set to the next value (currently beyond length of list)
        self.counter = self.size()

    def report_all(self, logger_name):
        assert logger_name is not None
        l = lg.getLogger(logger_name)
        # Get total number of warnings
        n = self.size()
        # No warnings? Great!
        if n == 0:
            l.subheading('Reported warnings')
            l('> No warnings!')
            return
        # Report warnings
        banner = '{} warnings/non-fatal errors (see below)'.format(n)
        l.subheading(banner)
        for i, e in enumerate(self.list()):
            l.bar()
            l('Warnings {} of {}'.format(i+1, n))
            l.bar()
            l(e)
        l.subheading(banner.replace('below','above'))


def get_warning_handler(logger, handler_name="warnings", recurse_parents=True):

    return get_handler_recursive(
        logger = logger,
        handler_name = handler_name,
        recurse_parents = recurse_parents,
        )

def get_handler_recursive(logger, handler_name, recurse_parents=True):
    """
    Looks for warning handlers in logger, and optionally parent loggers,
    until one is found or root is reached
    """

    handler_names = [h.name for h in logger.handlers]
    if handler_name in handler_names:
        index = handler_names.index(handler_name)
        return logger.handlers[index]

    # None found and is root: return none
    if (logger is logger.root) or (logger.name == 'root'):
        return None

    # Check parents of logger
    if (recurse_parents is True):
        return get_handler_recursive(
            logger = logger.parent,
            handler_name = handler_name,
            recurse_parents = recurse_parents,
        )

    return None

def setup_root_logging(formatter=None, level=lg.INFO):

    if formatter is None:
        formatter = lg.Formatter(fmt='%(message)s')

    # Get root logger
    logger = lg.getLogger()

    # Set logger to info by default
    logger.setLevel(level)

    # Add standard output to root logger
    if not logger.handlers:
        ch = lg.StreamHandler(stream=sys.stdout)
        ch.setFormatter(formatter)
        ch.setLevel(lg.DEBUG)
        logger.addHandler(ch)

    return logger

def get_logger_with_headings_maybe(name):
    """
    Returns LoggerWithHeadings named `name` if it hasn't been created,
    otherwise returns existing logger named `name`.
    """
    # Get the current class
    prevClass = lg.getLoggerClass()
    # Set the desired class
    lg.setLoggerClass(LoggerWithHeadings)
    # Get named logger
    logger = lg.getLogger(name)
    # Return logger class to previous
    lg.setLoggerClass(prevClass)

    return logger

# Shortcut for compatibility with python logging module
getLogger = get_logger_with_headings_maybe

def setup_logging_basic(name):
    """One liner to setup stdout logging (if not already done) and return a functioning Logger"""

    setup_root_logging()

    return get_logger_with_headings_maybe(name)

def add_warning_handler(logger, name='warnings'):

    # Custom warning collector
    wh = WarningListHandler(name=name)
    # Special warning formatter
    wfmt = lg.Formatter(fmt='%(levelname)s -- %(message)s')
    wh.setFormatter(wfmt)
    wh.setLevel(lg.WARNING)
    logger.addHandler(wh)

    return wh

def setup_logging(name, log_file=None, warning_handler_name='warnings', debug=False):
    """
    One liner to setup logging for a named scope with log file and warnings tracker.
    if name == '__main__' will set up the root logger.
    """

    setup_root_logging()

    # Plain formatter
    fmt = lg.Formatter(fmt='%(message)s')

    # Set level to at least info
    if debug is True:
        level = lg.DEBUG
    else:
        level = lg.INFO

    # Get logger (special case for name == '__main__')
    if name == '__main__':
        # Get the root logger
        logger = lg.getLogger()
    else:
        logger = get_logger_with_headings_maybe(name)

    # Set level of logger
    logger.setLevel(level)

    # Set output log file
    if (log_file is not None):
        fh = lg.FileHandler(log_file)
        fh.setFormatter(fmt)
        fh.setLevel(lg.DEBUG)
        logger.addHandler(fh)

    if (warning_handler_name is not None):

        wh = add_warning_handler(
            logger = logger, 
            name = warning_handler_name,
            )

    return get_logger_with_headings_maybe(name)

import sys
import logging as lg

def get_handler_recursive(logger, handler_name='warnings', recurse_parents=True):
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


class Bar:

    width = 40
    body = '-'
    head = '>>>'

    def __init__(self):
        self.body_width = max(0, self.width - len(self.head))
        self._bar = self.body*self.body_width + self.head

    def __call__(self, blank_before=False, blank_after=False):
        return '\n'*blank_before + self._bar + '\n'*blank_after


class Heading:

    width = 100
    spacer = '#'
    decorator = ' <~~~> '
    side_width = 3

    def __init__(self):
        self.content_width = max(0, self.width - (2*self.side_width))
        self.side_padding = self.side_width * self.spacer

    def __call__(self, text, spacer=False, blank=False):
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


class PandemicLogger(lg.Logger):

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

    def list(self):
        return iter(self.log_list)


class WarningListHandler(ListHandler):

    def __init__(self, name='warnings'):
        # run the regular ListHandler __init__
        ListHandler.__init__(self)
        # Custom attributes
        self.set_name(name)
        self.counter = 0

    def flush(self, logger_name=None):
        logger = lg.getLogger(logger_name)
        pass

    def report(self, logger_name=None):
        logger = lg.getLogger(logger_name)
        pass


def setup_root_logging(formatter=None, level=lg.INFO):

    if formatter is None:
        formatter = lg.Formatter(fmt='%(message)s')

    # Get root logger
    logger = lg.getLogger()

    # Add standard output to root logger
    ch = lg.StreamHandler(stream=sys.stdout)
    ch.setFormatter(formatter)
    ch.setLevel(level)
    logger.addHandler(ch)

    # Override base class to allow headings
    lg.setLoggerClass(PandemicLogger)

    return logger

def setup_logging(name, log_file=None, debug=False):
    """
    Setup logging for a named scope.
    if name == '__main__' will set up the root logger.
    """

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
        # Get a named logger
        logger = lg.getLogger(name)

    # Set level of logger
    logger.setLevel(level)

    # Set output log file
    if (log_file is not None):
        fh = lg.FileHandler(log_file)
        fh.setFormatter(fmt)
        fh.setLevel(level)
        logger.addHandler(fh)

    # Custom warning collector
    wh = WarningListHandler()
    # Special warning formatter
    wfmt = lg.Formatter(fmt='%(levelname)s -- %(message)s')
    wh.setFormatter(wfmt)
    wh.setLevel(lg.WARNING)
    logger.addHandler(wh)

    return logger

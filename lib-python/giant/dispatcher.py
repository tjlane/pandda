import giant.logs as lg
logger = lg.getLogger(__name__)

import procrunner
from giant.paths import is_available


class Dispatcher(object):

    def __init__(self, program):

        if not is_available(program):
            raise ValueError("Can't find the program '{!s}'. Is it available?".format(program))

        self.program = program
        self.command = [program]
        self.stdin = None
        self.working_directory = None
        self.timeout = None
        self.debug = False

        self.result = None

    def __str__(self):
        return (
            'Program:\n\t' + self.command[0] + '\n' +
            'Command line arguments:\n\t' + '\n\t'.join(self.command[1:] if self.command[1:] else ['None']) + '\n' +
            'Standard Input:\n\t' + (self.stdin if self.stdin is not None else 'None\n').replace('\n','\n\t')
        )

    def append_arg(self, arg):
        self.command.append(arg)

    def extend_args(self, args):
        self.command.extend(args)

    def append_stdin(self, line):
        if self.stdin is None:
            self.stdin = ''
        self.stdin += line.strip('\n') + '\n'

    def extend_stdin(self, lines):
        for line in lines:
            self.append_stdin(line)

    def as_string(self):
        return str(self)

    def as_command(self):
        out_str = ' '.join(self.command)
        if self.stdin:
            out_str += (' <<eof\n' + self.stdin + '\neof')
        return out_str

    def run(self):

        assert self.result is None, 'Program has already been run.'

        self.result = procrunner.run(
            command = self.command,
            timeout = self.timeout,
            debug = self.debug,
            stdin = (
                bytes(self.stdin, 'utf-8')
                if self.stdin is not None
                else None
                ),
            print_stdout = False,
            print_stderr = False,
            #callback_stdout = None,
            #callback_stderr = None,
            #environment = None,
            #environment_override = None,
            #win32resolve = True,
            working_directory = self.working_directory,
        )

        return self.result

    def write_output(self, log_file):

        separator = '\n----------------------\n'

        with open(log_file, 'a') as fh:

            fh.write(separator)
            fh.write('Command information:\n')
            fh.write(self.as_string()+'\n')
            fh.write(separator)
            fh.write('Command to re-run:\n')
            fh.write(self.as_command()+'\n')
            fh.write(separator)
            fh.write('Program STDOUT:\n')
            fh.write(str(self.result.stdout))
            fh.write(separator)
            fh.write('Program STDERR:\n')
            fh.write(str(self.result.stderr))
            fh.write(separator)



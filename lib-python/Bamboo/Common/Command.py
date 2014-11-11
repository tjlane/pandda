#! /usr/local/python/python2.7.3-64bit/bin/python

# Contains Functions for Launching Commands

import os, sys, subprocess, threading, shlex, time

# Code here adapted from code by Sebastian Kelm

class CommandManager(object):
    '''
        Enables commands be run on another thread, with various options

        Based on jcollado's solution:
        http://stackoverflow.com/questions/1191374/subprocess-with-timeout/4825933#4825933

        https://gist.github.com/1306188

        Modified by Sebastian Kelm: added timedout parameter
        Modified by Nicholas Pearce: added code inspired by code from the CCP4 dispatcher project
    '''

    def __init__(self, program):
        # Name of program
        self.program = shlex.split(program)
        self.command = self.program
        # Command line arguments
        self.args = []
        self.kwargs = {}
        # Process Control
        self.pipe = True
        self.timeout = None
        self.process = None
        self.timedout = False
        # Response
        self.inp = ''
        self.out = ''
        self.err = ''
        # Meta
        self.runtime = None

    def SetParameters(self, timeout=-1, pipe=-1):
        """Set various program parameters"""

        if timeout != -1: self.timeout = timeout
        if pipe != -1:    self.pipe = pipe

    def SetArguments(self, *args):
        """Set the CMD LINE inputs for the program - WILL OVERWRITE CURRENT ARGS"""

        self.args=[]
        for arg in args:
            # Takes any arguments given as a list 'as is', but splits any string arguments
            # Useful for passing strings to it => pass string inside a list
            if isinstance(arg,list):
                self.args.extend(arg)
            elif isinstance(arg,str):
                self.args.extend(shlex.split(arg))
        # Add to the Command
        self.command = self.program + self.args

    def SetInput(self, input):
        """Set the STDIN for the program - WILL OVERWRITE CURRENT INPUT"""

        # Format the input correctly
        if isinstance(input,list):
            input = "\n".join(input)
        elif not isinstance(input, str):
            raise TypeError('input must be str or list')
        if not input.endswith("\n"):
            input = input + "\n"
        self.inp = input

    def AppendInput(self, input):
        """Append the STDIN for the program"""

        # Format
        if isinstance(input,list):
            input = "\n".join(input)
        if not input.endswith("\n"):
            input = input + "\n"
        self.inp += input

    def PrintSettings(self):
        """Print out the current settings of the object"""

        print 'Command:', self.command
        print 'Inputs:', self.inp

    def _PrepareInputs(self):
        """Prepare the Input Pipes"""

        if self.inp.strip():
            self.__stdin__ = subprocess.PIPE
            self.kwargs['stdin'] = self.__stdin__
        else:
            self.__stdin__ = None
            self.kwargs.pop('stdin',None)

    def _PrepareOutputs(self):
        """Prepare the Output Pipes"""

        if self.pipe:
            self.__stdout__ = subprocess.PIPE
            self.__stderr__ = subprocess.PIPE
            self.kwargs['stdout'] = self.__stdout__
            self.kwargs['stderr'] = self.__stderr__

    def Run(self):
        """Run command with keyword arguments"""

        self._PrepareInputs()
        self._PrepareOutputs()

        # Function to run the subprocess
        def target():
            self.process = subprocess.Popen(args=self.command, **self.kwargs)
            self.out, self.err = self.process.communicate(self.inp)

        # Create the new thread
        self.thread = threading.Thread(target=target)
        # Start time
        self.t1 = time.time()
        # Start the Thread
        self.thread.start()
        # Check for timeout
        self.thread.join(self.timeout)
        # Act accordingly
        if self.thread.is_alive():
            self.process.terminate()
            self.thread.join()
            self.timedout = True
            return 666
        # End time
        self.t2 = time.time()
        # Calculate runtime
        self.runtime = self.t2-self.t1
        # Return process exit status
        return self.process.returncode

    def run_async(self):
        """Run Command with the Keyword Arguments and then return the thread object"""

        self._PrepareInputs()
        self._PrepareOutputs()

        # Function to run the subprocess
        def target():
            self.process = subprocess.Popen(args=self.command, **self.kwargs)
            self.out, self.err = self.process.communicate(self.inp)

        # Create the new thread
        self.thread = threading.Thread(target=target)
        # Start time
        self.t1 = time.time()
        # Start the Thread
        self.thread.start()

    def join_async(self):
        """Join the thread object returned by self.run_async"""

        # Check for timeout
        self.thread.join(self.timeout)
        # Act accordingly
        if self.thread.is_alive():
            self.process.terminate()
            self.thread.join()
            self.timedout = True
            return 666
        # End time
        self.t2 = time.time()
        # Calculate runtime
        self.runtime = self.t2-self.t1
        # Return process exit status
        return self.process.returncode

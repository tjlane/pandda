#! /usr/local/python/python2.7.3-64bit/bin/python

from argparse import ArgumentParser

class argParser(ArgumentParser):
    """Minor Modifications to Argument Parser"""
    def parseArguments(self, args=None):
        """Parse Arguments and return results also as a dictionary"""
        if args:
            self.Inputs = self.parse_args(args)
        else:
            self.Inputs = self.parse_args()

        return self.Inputs, self.Inputs.__dict__


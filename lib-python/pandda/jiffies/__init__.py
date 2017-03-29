import os, sys, copy

from pandda import HEADER_TEXT

# This will fail inside coot...
try:
    from giant.jiffies import run_default
    # Override the default values of the run_default class
    run_default.HEADER_TEXT = HEADER_TEXT
except:
    pass


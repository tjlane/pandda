import os, sys, copy

from giant.jiffies import run_default

from pandda import HEADER_TEXT

# Override the default values of the run_default class
run_default.HEADER_TEXT = HEADER_TEXT

import os, sys, copy

from pandda import module_info

# This will fail inside coot...
try:
    from giant.jiffies import run_default
    # Override the default values of the run_default class
    run_default._module_info = module_info
except:
    pass


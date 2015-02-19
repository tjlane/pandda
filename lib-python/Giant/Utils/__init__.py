

import os, sys, glob

def status_bar(n, n_max):
    n_jump = max(1,int(n_max/100))
    if (n+1==n_max):
        print '\r>>', 100, '%'
    elif (n%n_jump==0):
        print '\r>>', int(round(100.0*(n+1)/n_max,0)), '%',
    sys.stdout.flush()


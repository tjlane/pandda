

import os, sys, glob

def status_bar(n, n_max):
    n_jump = max(1,int(n_max/100))
    assert n<=n_max, 'n must be less than or equal to n_max'
    if (n==n_max):
        print '\r>>', 100, '%'
    elif (n%n_jump==0):
        print '\r>>', int(round(100.0*n/n_max,0)), '%',
    sys.stdout.flush()


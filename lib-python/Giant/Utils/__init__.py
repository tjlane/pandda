

import os, sys, glob

def status_bar(n, n_max):
    n_jump = max(1,int(n_max/100))
    assert n<=n_max, 'n must be less than or equal to n_max'
    if (n==n_max):
        print '\r>>', 100, '%'
    elif (n%n_jump==0):
        print '\r>>', int(round(100.0*n/n_max,0)), '%',
    sys.stdout.flush()

def rel_symlink(orig, link):
    """Make a relative symlink from link to orig"""
    assert os.path.exists(orig), 'FILE DOES NOT EXIST: {!s}'.format(orig)
    assert not os.path.exists(link), 'LINK ALREADY EXISTS: {!s}'.format(link)
    orig = os.path.abspath(orig)
    link = os.path.abspath(link)
    assert not link.endswith('/'), 'LINK CANNOT END WITH /'
    os.symlink(os.path.relpath(orig, start=os.path.dirname(link)), link)



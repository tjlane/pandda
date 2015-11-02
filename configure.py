import os, sys

from distutils.spawn import find_executable

##########################################################
# Find the path of the running python
##########################################################
this_python = sys.executable

##########################################################
# Get the path of this script
##########################################################
script_path = os.path.realpath(__file__)
print 'RUNNING: {}'.format(script_path)
##########################################################
# Get the path of the PANDDAs top folder
##########################################################
PANDDA_TOP = os.path.dirname(script_path)
print 'PANDDAs RUNNING FROM: {}'.format(PANDDA_TOP)
##########################################################
# Find the bin of the PANDDAs module
##########################################################
PANDDA_BIN = os.path.join(PANDDA_TOP, 'bin')
assert os.path.exists(PANDDA_BIN)

##########################################################
# Get the user-given path to the cctbx.python to be used
##########################################################
cctbx_python = [arg for arg in sys.argv if 'cctbx.python' in arg]
pandda_python = os.path.join(PANDDA_BIN, 'pandda.python')
if not cctbx_python:
    raise Exception('No cctbx path given - pass path to cctbx.python executable')
else:
    cctbx_python = find_executable(cctbx_python[0])
# Link cctbx.python to pandda.python
if not cctbx_python:
    raise Exception('No cctbx path found - pass path to cctbx.python executable')
else:
    print('\n=====================================+>')
    print('Updating pandda.python paths')
    print('=====================================+>')
    # Create a link to this python
    if os.path.islink(pandda_python):
        print 'Unlinking: {}'.format(pandda_python)
        os.unlink(pandda_python)
    elif os.path.exists(pandda_python):
        raise Exception('{} is hardcoded (should be a symbolic link) - cannot replace'.format(pandda_python))
    # Link the python to a new pandda.python symlink
    print 'LINKING {} -> {}'.format(cctbx_python, pandda_python)
    os.symlink(cctbx_python, pandda_python)

##########################################################
# Get the user-given path to the libtbx.ipython to be used
##########################################################
cctbx_ipython = [arg for arg in sys.argv if 'libtbx.ipython' in arg]
pandda_ipython = os.path.join(PANDDA_BIN, 'pandda.ipython')
if not cctbx_ipython:
    print('No libtbx.ipython path given - pass path to libtbx.ipython executable')
else:
    cctbx_ipython = find_executable(cctbx_ipython[0])
# Link libtbx.ipython to pandda.ipython
if not cctbx_ipython:
    print('No libtbx.ipython path found - pass path to libtbx.ipython executable')
else:
    print('\n=====================================+>')
    print('Updating pandda.ipython paths')
    print('=====================================+>')
    # Create a link to this ipython
    if os.path.islink(pandda_ipython):
        print 'Unlinking: {}'.format(pandda_ipython)
        os.unlink(pandda_ipython)
    elif os.path.exists(pandda_ipython):
        raise Exception('{} is hardcoded (should be a symbolic link) - cannot replace'.format(pandda_ipython))
    # Link the ipython to a new pandda.ipython symlink
    print 'LINKING {} -> {}'.format(cctbx_ipython, pandda_ipython)
    os.symlink(cctbx_ipython, pandda_ipython)

##########################################################
# Print lines to be added to bashrc, cshrc
##########################################################
print('\n=====================================+>')
print('Add following lines to bashrc:')
print('=====================================+>')
print('export PANDDA_TOP={}'.format(PANDDA_TOP))
print('source {}'.format(os.path.join('$PANDDA_TOP','scripts','pandda.sh')))
print('\n=====================================+>')
print('Add following lines to cshrc:')
print('=====================================+>')
print('setenv PANDDA_TOP {}'.format(PANDDA_TOP))
print('source {}'.format(os.path.join('$PANDDA_TOP','scripts','pandda.csh')))
print('\n=====================================+>')





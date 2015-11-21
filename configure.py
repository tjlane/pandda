import os, sys

from distutils.spawn import find_executable

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

PANDDA_BIN = os.path.join(PANDDA_TOP, 'bin')
assert os.path.exists(PANDDA_BIN)

pandda_python  = os.path.join(PANDDA_BIN, 'pandda.python')
pandda_ipython = os.path.join(PANDDA_BIN, 'pandda.ipython')
pandda_sh      = os.path.join(PANDDA_TOP, 'scripts', 'pandda.sh')
pandda_csh     = os.path.join(PANDDA_TOP, 'scripts', 'pandda.csh')

##########################################################
# Get the user-given path to the cctbx.python to be used
##########################################################
cctbx_python = [arg for arg in sys.argv if 'cctbx.python' in arg]
if cctbx_python: cctbx_python = find_executable(cctbx_python[0])
else:            cctbx_python = find_executable('cctbx.python')

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
if cctbx_ipython: cctbx_ipython = find_executable(cctbx_ipython[0])
else:             cctbx_ipython = find_executable('libtbx.ipython')

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
# Write bashrc
##########################################################
with open(pandda_sh, 'w') as fh:
    fh.write('# Setup for PANDDAs\n')
    fh.write('export PANDDA_TOP={}\n'.format(PANDDA_TOP))
    fh.write('export PYTHONPATH=${PYTHONPATH}:$PANDDA_TOP/lib-python\n')
    fh.write('export PATH=${PATH}:$PANDDA_TOP/bin\n')

##########################################################
# Write cshrc
##########################################################
with open(pandda_csh, 'w') as fh:
    fh.write('# Setup for PANDDAs\n')
    fh.write('setenv PANDDA_TOP {}\n'.format(PANDDA_TOP))
    fh.write('setenv PYTHONPATH ${PYTHONPATH}:$PANDDA_TOP/lib-python\n')
    fh.write('setenv PATH ${PATH}:$PANDDA_TOP/bin\n')

##########################################################
# Print lines to be added to bashrc, cshrc
##########################################################
print('\n=====================================+>')
print('Add following lines to bashrc:')
print('=====================================+>')
print('source {}'.format(pandda_sh))
print('\n=====================================+>')
print('Add following lines to cshrc:')
print('=====================================+>')
print('source {}'.format(pandda_csh))
print('\n=====================================+>')





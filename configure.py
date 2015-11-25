import os, sys

from distutils.spawn import find_executable

##########################################################
# Get the path of this script
##########################################################
script_path = os.path.realpath(__file__)
print('=====================================+>')
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
pandda_coot    = os.path.join(PANDDA_BIN, 'pandda.coot')

pandda_sh      = os.path.join(PANDDA_TOP, 'scripts', 'pandda.sh')
pandda_csh     = os.path.join(PANDDA_TOP, 'scripts', 'pandda.csh')

##########################################################
# Get the user-given path to the cctbx.python to be used
##########################################################
cctbx_python = [arg for arg in sys.argv if 'cctbx.python' in arg]
if cctbx_python: cctbx_python = find_executable(cctbx_python[0])
else:            cctbx_python = None

if not cctbx_python:
    print('=====================================+>')
    if os.path.exists(pandda_python):
        print('Not updating python link - pandda.python already exists')
        print('Pass argument `cctbx.python` to override')
        cctbx_python = None
    else:
        print('Will try to use default cctbx.python')
        cctbx_python = find_executable('cctbx.python')

if cctbx_python:
    print('=====================================+>')
    print('Updating pandda.python paths:')
    # Create a link to this python
    if os.path.islink(pandda_python):
        print 'Unlinking: {}'.format(pandda_python)
        os.unlink(pandda_python)
    elif os.path.exists(pandda_python):
        raise Exception('{} is hardcoded (should be a symbolic link) - cannot replace'.format(pandda_python))
    # Link the python to a new pandda.python symlink
    print 'LINKING {} -> {}'.format(cctbx_python, pandda_python)
    os.symlink(cctbx_python, pandda_python)
elif os.path.exists(pandda_python):
    pass
else:
    print('=====================================+>')
    raise Exception('No cctbx path found - pass path to cctbx.python executable')

##########################################################
# Get the user-given path to the libtbx.ipython to be used
##########################################################
cctbx_ipython = [arg for arg in sys.argv if 'libtbx.ipython' in arg]
if cctbx_ipython: cctbx_ipython = find_executable(cctbx_ipython[0])
else:             cctbx_ipython = None

# Link libtbx.ipython to pandda.ipython
if not cctbx_ipython:
    print('=====================================+>')
    if os.path.exists(pandda_ipython):
        print('Not updating ipython link - pandda.ipython already exists')
        print('Pass argument `libtbx.ipython` to override')
        cctbx_ipython = None
    else:
        print('Will try to use default libtbx.ipython')
        cctbx_ipython = find_executable('libtbx.ipython')

if cctbx_ipython:
    print('=====================================+>')
    print('Updating pandda.ipython paths:')
    # Create a link to this ipython
    if os.path.islink(pandda_ipython):
        print 'Unlinking: {}'.format(pandda_ipython)
        os.unlink(pandda_ipython)
    elif os.path.exists(pandda_ipython):
        raise Exception('{} is hardcoded (should be a symbolic link) - cannot replace'.format(pandda_ipython))
    # Link the ipython to a new pandda.ipython symlink
    print 'LINKING {} -> {}'.format(cctbx_ipython, pandda_ipython)
    os.symlink(cctbx_ipython, pandda_ipython)
elif os.path.exists(pandda_ipython):
    pass
else:
    print('=====================================+>')
    print('No libtbx.ipython path found - pass path to libtbx.ipython executable to link')

##########################################################
# Create pandda.coot script
##########################################################
ld_path = [arg for arg in sys.argv if arg.endswith('lib')]
if ld_path: ld_path = os.path.realpath(ld_path[0])
else:       ld_path = None
coot_path = [arg for arg in sys.argv if 'coot' in arg]
if coot_path: coot_path = find_executable(coot_path[0])
else:         coot_path = None

if not coot_path:
    if os.path.exists(pandda_coot):
        print('=====================================+>')
        print('Not updating coot - pandda.coot already exists')
        print('Pass argument `coot` to override')
        coot_path = None
    else:
        print('=====================================+>')
        print('Will try to use default coot')
        coot_path = find_executable('coot')

if coot_path:
    print('=====================================+>')
    print('Creating pandda.coot script:')
    print('Using coot from {}'.format(coot_path))
    if ld_path: print('Using LD_PATH of {}'.format(ld_path))
    with open(pandda_coot, 'w') as fh:
        fh.write('#!/usr/bin/env sh\n')
        if ld_path: fh.write('LD_LIBRARY_PATH={} '.format(ld_path))
        fh.write('{} --python $@\n'.format(coot_path))
    os.chmod(pandda_coot, 0751)
elif os.path.exists(pandda_coot):
    pass
else:
    print('=====================================+>')
    print('No coot path found - is coot installed or findable in path?')

##########################################################
print('=====================================+>')
##########################################################
# Write bashrc
##########################################################
print('Writing pandda bash script.')
with open(pandda_sh, 'w') as fh:
    fh.write('# Setup for PANDDAs\n')
    fh.write('export PANDDA_TOP={}\n'.format(PANDDA_TOP))
    fh.write('export PYTHONPATH=${PYTHONPATH}:$PANDDA_TOP/lib-python\n')
    fh.write('export PATH=${PATH}:$PANDDA_TOP/bin\n')
##########################################################
# Write cshrc
##########################################################
print('Writing pandda csh script.')
with open(pandda_csh, 'w') as fh:
    fh.write('# Setup for PANDDAs\n')
    fh.write('setenv PANDDA_TOP {}\n'.format(PANDDA_TOP))
    fh.write('setenv PYTHONPATH ${PYTHONPATH}:$PANDDA_TOP/lib-python\n')
    fh.write('setenv PATH ${PATH}:$PANDDA_TOP/bin\n')

##########################################################
# Print lines to be added to bashrc, cshrc
##########################################################

print('')
print('=====================================+>')
print('Add following lines to bashrc:')
print('=====================================+>')
print('source {}'.format(pandda_sh))
print('=====================================+>')
print('Add following lines to cshrc:')
print('=====================================+>')
print('source {}'.format(pandda_csh))
print('=====================================+>')





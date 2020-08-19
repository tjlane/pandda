import os, sys, glob
from setuptools import setup, find_packages, findall
from distutils.spawn import find_executable

import imp
giant = imp.load_source('giant', os.path.join(os.path.dirname(__file__), 'lib-python/giant/__init__.py'))
VERSION = giant.__version__

#####################################################################################

def requirements():
    # The dependencies are the same as the contents of requirements.txt
    with open('requirements.txt') as f:
        return [line.strip() for line in f if line.strip()]

#####################################################################################

if '--for-ccp4' in sys.argv:
    # Select scripts for install
    install_scripts = [ 'bin/pandda.analyse',
                        'bin/pandda.inspect',
                        'bin/pandda.export',
                        'bin/giant.score_model',
                        'bin/giant.quick_refine',
                        'bin/giant.make_restraints',
                        'bin/giant.merge_conformations',
                        'bin/giant.split_conformations',
                      ]
else:
    # Standard install
    install_scripts = findall(dir='bin')

# Enable the python command to be changed
if '--python' in sys.argv:
    python = sys.argv[sys.argv.index('--python')+1]
    # And clean up
    sys.argv.remove('--python')
    sys.argv.remove(python)
elif '--for-ccp4' in sys.argv:
    python = 'ccp4-python'
else:
    python = 'cctbx.python'

# And clean up
try: sys.argv.remove('--for-ccp4')
except: pass

assert find_executable(python), "Can't find executable - is this the correct python: {}".format(python)

# Modify the scripts to use the correct python/coot
for s in install_scripts:
    with open(s, 'r') as fh: s_conts = fh.read()
    s_conts = s_conts.replace('cctbx.python', python)
    with open(s, 'w') as fh: fh.write(s_conts)

#####################################################################################

setup(
    name = 'panddas',
    version = VERSION,
    description = 'Multi-dataset crystallographic analyses',
    author = 'Nicholas M Pearce',
    author_email = 'nicholas.pearce.0@gmail.com',
    url = 'http://pandda.bitbucket.org',
    license = 'CC BY-SA 4.0',
    install_requires = requirements(),
    package_dir = {'':'lib-python'},
    scripts = install_scripts,
    packages = find_packages(
        where = 'lib-python',
        exclude = (
            'prototypes*',
            '*egg*',
            ),
        ),
    include_package_data= True,
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Operating System :: Unix',
        ],
)

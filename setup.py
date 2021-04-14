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

# Standard install
install_scripts = findall(dir='bin')

#####################################################################################

setup(
    name = 'panddas',
    version = VERSION,
    description = 'Multi-dataset crystallographic analyses',
    author = 'Nicholas M Pearce',
    author_email = 'nicholas.pearce.0@gmail.com',
    url = 'http://pandda.bitbucket.org',
    license = 'LGPL-3.0-or-later',
    #license_files = ('LICENSE.txt',),
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

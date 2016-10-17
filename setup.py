import os, sys, glob
from setuptools import setup, find_packages, findall
from distutils.spawn import find_executable

#####################################################################################

if '--for-ccp4' in sys.argv:
    # Select scripts for install
    install_scripts = [ 'bin/pandda.analyse',
                        'bin/pandda.inspect',
                        'bin/pandda.export',
                        'bin/giant.score_model',
                        'bin/giant.quick_refine',
                        'bin/giant.create_occupancy_params',
                        'bin/giant.strip_conformations'         ]
    # And clean up
    sys.argv.remove('--for-ccp4')
else:
    # Standard install
    install_scripts = findall(dir='bin')

# Enable the python command to be changed
if '--python' in sys.argv:
    python = sys.argv[sys.argv.index('--python')+1]
    # And clean up
    sys.argv.remove('--python')
    sys.argv.remove(python)
else:
    python = 'ccp4-python'

print('PanDDA scripts will run using: {}'.format(python))
assert find_executable(python), "Can't find executable - is this the correct python: {}".format(python)

# Modify the scripts to use the correct python/coot
for s in install_scripts:
    with open(s, 'r') as fh: s_conts = fh.read()
    s_conts = s_conts.replace('pandda.python', python)
    s_conts = s_conts.replace('pandda.coot',   'coot')
    with open(s, 'w') as fh: fh.write(s_conts)

#####################################################################################

setup(  name                = 'panddas',
        version             = '0.1.11',
        description         = 'Multi-dataset crystallographic analyses',
        author              = 'Nicholas M Pearce',
        author_email        = 'nicholas.pearce.0@gmail.com',
        url                 = 'http://pandda.bitbucket.org',
        license             = 'CC BY-SA 4.0',
        install_requires    = ( 'pandas',
                                'jinja2',
                                'markupsafe',
                                'ascii_graph',
                                ),
        package_dir         = {'':'lib-python'},
        scripts             = install_scripts,
        packages            = find_packages(where   ='lib-python',
                                            exclude = ( 'pandemic',     'pandemic.*',
                                                        'phenix_pandda*',
                                                        'prototypes*',
                                                        '*egg*',
                                                        )
                                            ),
        include_package_data= True,
        classifiers         = ( 'Development Status :: 4 - Beta',
                                'Intended Audience :: Developers',
                                'Natural Language :: English',
                                'Programming Language :: Python',
                                'Programming Language :: Python :: 2.7',
                                'Topic :: Software Development :: Libraries :: Python Modules',
                                'Operating System :: Unix',
                                ),
        )

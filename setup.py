import os, sys, glob
from setuptools import setup, find_packages, findall

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
    # Modify the scripts for ccp4
    for s in install_scripts:
        with open(s, 'r') as fh: s_conts = fh.read()
        s_conts = s_conts.replace('pandda.python', 'ccp4-python')
        s_conts = s_conts.replace('pandda.coot',   'coot')
        with open(s, 'w') as fh: fh.write(s_conts)
    # And clean up
    sys.argv.remove('--for-ccp4')
else:
    # Standard install
    install_scripts = findall(dir='bin')

#####################################################################################

setup(  name                = 'panddas',
        version             = '0.1.9',
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

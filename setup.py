import os, sys, glob
from setuptools import setup, find_packages, findall

setup(  name                = 'panddas',
        version             = '1.0.0',
        description         = 'Multi-dataset crystallographic analyses',
        author              = 'Nicholas M Pearce',
        author_email        = 'nicholas.pearce.0@gmail.com',
        url                 = 'pandda.bitbucket.org',
        license             = 'CC BY-SA 4.0',
        install_requires    = ( 'pandas',
                                'jinja2',
                                'markupsafe',
                                'ascii_graph',
                                ),
        package_dir         = {'':'lib-python'},
        scripts             = findall(dir='bin'),
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

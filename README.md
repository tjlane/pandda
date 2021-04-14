
## System Requirements

1. This package is written for Mac/Linux. It is not guaranteed to be compatible with Windows.
 
2. This package must be installed within an up-to-date ccp4/phenix/cctbx environment (cctbx version v2021.1 or later)
    * (recommended) For Phenix, the latest version can be found here: [Phenix Download Page](https://www.phenix-online.org/download/)
    * For CCP4, the latest version can be found here: [CCP4 Download Page](http://www.ccp4.ac.uk/download/)

## Setup instructions

Typical install time: `~1 minute`.

There are currently some difficulties and extra required steps to install the panddas package mostly due to the demise and somewhat zombification of python2.7. Hopefully things will be simpler in the world of python3.
However, for the moment: 

### Install in PHENIX (RECOMMENDED)

1. To install the newest version of panddas, run
    * ```git clone https://bitbucket.org/pandda/pandda.git ; cd pandda ; cctbx.python setup.py install```

2. The scripts will be installed into a directory that is not in the standard PHENIX PATH variables
    * Therefore you will also need to run
        * BASH/SH users: ```export PATH=$PHENIX/conda_base/bin:$PATH```
        * CSH users: ```setenv PATH $PHENIX/conda_base/bin:$PATH```
    * You may want to add the appropriate line above to your rc file (e.g. `~/.bashrc`).
  
3. To check that everything has installed properly, you should be able to run e.g. ```pandemic.adp "?"``` and be given a list of options. To check the functionality of each of the programs, follow one of the tutorials on the webpage in the section below.

### Install in CCP4

1. Check that pip is installed
    * Run
        * ```ccp4-python -m pip --version```
    * If you get a message "No module named pip", then pip is not installed, and you should run 
        * ```curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py ; ccp4-python get-pip.py```

2. Uninstall any existing versions of panddas (as they will not always uninstall automatically) 
    * ```ccp4-python -m pip uninstall panddas -y``` 

3. To install the newest version of panddas, run
    * ```git clone https://bitbucket.org/pandda/pandda.git ; cd pandda ; ccp4-python setup.py install```

4. Depending on your system, this may install the scripts folder that is not added to the PATH by the CCP4 setup script. 
    * To identify if this is the case, run 
        * ```which pandemic.adp```
    * If this does not find the script in the new CCP4 install then you need to add a new folder to your path.
        * First, run 
            * ```find $CCP4 -name pandemic.adp```.
            * This should reveal two scripts, one of which is in a folder named "bin".
        * Next, add the bin folder that contains the script to the PATH
            * on BASH/SH: ```export PATH=<path>:$PATH```
            * on CSH: ```setenv PATH <path>:$PATH```
            * where <path> is the folder containing the pandemic.adp script from above.
        * You may want to add the appropriate line above to your rc file (e.g. `~/.bashrc`).

5. To check that everything has installed properly, you should be able to run  e.g. ```pandemic.adp "?"``` and be given a list of options. To check the functionality of each of the programs, follow one of the tutorials on the webpage in the section below.

### Troubleshooting

* The panddas scripts use the first available cctbx.python available in the PATH to run the programs. Therefore, if you have another sourced CCP4/PHENIX distribution available, this may get used instead, causing errors. Therefore, make sure you are not sourcing another distribution when trying to use panddas. Alternatively, you may simply resource the CCP4/PHENIX setup script for the distribution that you want to use and this will fix most problems.

## Tutorials and Demos

* Tutorials and example output can be found of at [https://pandda.bitbucket.io](https://pandda.bitbucket.io)

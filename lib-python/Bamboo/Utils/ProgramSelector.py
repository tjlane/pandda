#! /usr/local/python/python2.7.3-64bit/bin/python

# Lists/Constants

# =================== #
# Refiner Settings    #
# =================== #

refiners = ['refmac','phenix.refine']
allowed_refiner_args = []

# =================== #
# Ligand Generation   #
# =================== #

ligand_builders = ['elbow','grade','writedict']
#ligand_builders = ['elbow','grade','prodrg','writedict']
allowed_builder_args = []

# =================== #
# Ligand Fitting      #
# =================== #

ligand_fitters = ['ligandfit','rhofit','flynn']
#ligand_fitters = ['phenix.ligandfit','arpwarp','rhofit','coot','jligand','flynn','afitt']
allowed_fitter_args = []

# =================== #

def get_ligand_builder(program):
    """Get Ligand Builder Object for `program`"""

    # Check it's a valid program
    if program not in ligand_builders:
        raise LigandBuilderSelectionError('{!s} not in {!s}'.format(program, ', '.join(ligand_builders)))

    if program in ['elbow','phenix.elbow']:
        from Bamboo.Wrappers.Builders.Elbow import elbowObject
        return elbowObject()
    elif program == 'grade':
        from Bamboo.Wrappers.Builders.Grade import gradeObject
        return gradeObject()
    elif program == 'prodrg':
        from Bamboo.Wrappers.Builders.Prodrg import prodrgObject
        return prodrgObject()
    elif program in ['writedict','afitt']:
        from Bamboo.Wrappers.Builders.Writedict import writedictObject
        return writedictObject()
    else:
        raise LigandBuilderSelectionError('MY BAD. The Code for {!s} has not been written yet...'.format(program))

def get_ligand_fitter(program):
    """Get Ligand Fitter Object for `program`"""

    # Check it's a valid program
    if program not in ligand_fitters:
        raise LigandFitterSelectionError('{!s} not in {!s}'.format(program, ', '.join(ligand_fitters)))

    if program in ['ligandfit','phenix.ligandfit']:
        from Bamboo.Wrappers.Fitters.LigandFit import ligandfitObject
        return ligandfitObject()
    elif program == 'arpwarp':
        from Bamboo.Wrappers.Fitters.ARPwARP import arpwarpObject
        return arpwarpObject()
    elif program == 'rhofit':
        from Bamboo.Wrappers.Fitters.Rhofit import rhofitObject
        return rhofitObject()
    elif program in ['flynn','afitt']:
        from Bamboo.Wrappers.Fitters.Flynn import flynnObject
        return flynnObject()
    elif program == 'coot':
        from Bamboo.Wrappers.Fitters.Coot import cootObject
        return cootObject()
    else:
        raise LigandFitterSelectionError('MY BAD. The Code for {!s} has not been written yet...'.format(program))

def get_refiner(program):
    '''Returns Initialised Refiner Object'''

    # Check it's a valid program
    if program not in refiners:
        raise RefinerSelectionError('{!s} not in {!s}'.format(program, ', '.join(refiners)))

    if program.lower()=='phenix.refine':
        # Get the phenix.refine refiner class
        from Bamboo.Wrappers.Refiners.Phenix import phenixrefineObject
        return phenixrefineObject()
    elif program.lower()=='refmac':
        # Get the refmac refiner class
        from Bamboo.Wrappers.Refiners.Refmac import refmacObject
        return refmacObject()
    else:
        # Not Found
        raise RefinerSelectionError('MY BAD. The Code for {!s} has not been written yet...'.format(program))


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
        from bamboo.wrappers.ligand_builders.elbow import ElbowObject
        return ElbowObject()
    elif program == 'grade':
        from bamboo.wrappers.ligand_builders.grade import GradeObject
        return GradeObject()
    elif program == 'prodrg':
        from bamboo.wrappers.ligand_builders.prodrg import ProdrgObject
        return ProdrgObject()
    elif program in ['writedict','afitt']:
        from bamboo.wrappers.ligand_builders.writedict import WritedictObject
        return WritedictObject()
    else:
        raise LigandBuilderSelectionError('MY BAD. The Code for {!s} has not been written yet...'.format(program))

def get_ligand_fitter(program):
    """Get Ligand Fitter Object for `program`"""

    # Check it's a valid program
    if program not in ligand_fitters:
        raise LigandFitterSelectionError('{!s} not in {!s}'.format(program, ', '.join(ligand_fitters)))

    if program in ['ligandfit','phenix.ligandfit']:
        from bamboo.wrappers.ligand_fitters.ligandfit import LigandfitObject
        return LigandfitObject()
    elif program == 'arpwarp':
        from bamboo.wrappers.ligand_fitters.arpwarp import ArpwarpObject
        return ArpwarpObject()
    elif program == 'rhofit':
        from bamboo.wrappers.ligand_fitters.rhofit import RhofitObject
        return RhofitObject()
    elif program in ['flynn','afitt']:
        from bamboo.wrappers.ligand_fitters.flynn import FlynnObject
        return FlynnObject()
    elif program == 'coot':
        from bamboo.wrappers.ligand_fitters.coot import CootObject
        return CootObject()
    else:
        raise LigandFitterSelectionError('MY BAD. The Code for {!s} has not been written yet...'.format(program))

def get_refiner(program):
    '''Returns Initialised Refiner Object'''

    # Check it's a valid program
    if program not in refiners:
        raise RefinerSelectionError('{!s} not in {!s}'.format(program, ', '.join(refiners)))

    if program.lower()=='phenix.refine':
        # Get the phenix.refine refiner class
        from bamboo.wrappers.refiners.phenix import PhenixrefineObject
        return PhenixrefineObject()
    elif program.lower()=='refmac':
        # Get the refmac refiner class
        from bamboo.wrappers.refiners.refmac import RefmacObject
        return RefmacObject()
    else:
        # Not Found
        raise RefinerSelectionError('MY BAD. The Code for {!s} has not been written yet...'.format(program))


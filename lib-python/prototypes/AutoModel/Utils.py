
from bamboo.common.parser import argParser
from bamboo.wrappers.program_selector import ligand_builders, ligand_fitters, refiners
from bamboo.common.logs import LigandRecord, PdbRecord, MtzRecord


# Ligand Model Record
class modelLog(object):
    """Object to store information about the ligand modelling process"""

    def __init__(self, ligsmile, refpdb, apopdb, inputmtz, builder, fitter, refiner, outdir, name='ModelLog', verbose=False, **kwargs):

        # Program used to generate the ligand model
        self.builder = builder
        # Program used to fit the ligand model to electron density
        self.fitter = fitter
        # Program used to refine the structures
        self.refiner = refiner
        # Directory of output files
        self.outdir = outdir

        # Record Objects
        self.lig = LigandRecord(name=name+'-Ligand', builder=builder, smile=ligsmile)
        self.mtz = MtzRecord(name=name+'-MTZ')
        if refpdb:
            self.pdb = PdbRecord(name=name+'-PDB-FROM REFERENCE')
            self.setReferencePDB(refpdb)
            self.setRawMTZ(inputmtz)
            self.pdb.history.setOriginalFile(refpdb)
        elif apopdb:
            self.pdb = PdbRecord(name=name+'-PDB-FROM APO FILE')
            self.setApoPDB(apopdb)
            self.setApoMTZ(inputmtz)
            self.pdb.history.setOriginalFile(apopdb)
        # Set MTZ to original
        self.mtz.history.setOriginalFile(inputmtz)

        # Set verboseness
        self.verbose = verbose

    def __str__(self):
        return 'ModelLog: {!s}'.format(self.name)

    def setReferencePDB(self, File, tag=None):
        if not self.pdb.reference:
            self.pdb.reference = File
            if tag: self.pdb.reference.tag=tag
        else:
            raise AssignmentError('Cannot Overwrite Reference PDB')
    def setRawMTZ(self, File, tag=None):
        if not self.mtz.raw:
            self.mtz.raw = File
            if tag: self.mtz.raw.tag=tag
        else:
            raise AssignmentError('Cannot Overwrite Raw MTZ')
    def setApoPDB(self, File, tag=None):
        if not self.pdb.apo:
            self.pdb.apo = File
            if tag: self.pdb.apo.tag=tag
        else:
            raise AssignmentError('Cannot Overwrite Apo PDB')
    def setApoMTZ(self, File, tag=None):
        if not self.mtz.apo:
            self.mtz.apo = File
            if tag: self.mtz.apo.tag=tag
        else:
            raise AssignmentError('Cannot Overwrite Apo MTZ')
    def addNewPDB(self, File, tag=None):
        self.pdb.current = File
        if tag: self.pdb.current.tag=tag
    def addNewMTZ(self, File, tag=None):
        self.mtz.current = File
        if tag: self.mtz.current.tag=tag

def modelLigandParser():
    """Parser for Ligand-Modelling Program"""

    des = """Parser for use with Ligand Modelling Pipelines"""
    # Initialise main parser object
    p = argParser(description=des)
    # Add arguments
    p.add_argument('--verbose','-v', action='store_true', default=False, help='Set Verboseness')
    # Ligand Arguements
    p.add_argument('--ligpdb', metavar='FILE', help='Ligand PDB file')
    p.add_argument('--ligcif', metavar='FILE', help='Ligand CIF file')
    p.add_argument('--ligsmile','--smile', help='Smile of ligand to be fitted')
    # Protein Arguments
    p.add_argument('--apopdb','--apo', metavar='FILE', help='Apo structure of the protein')
    p.add_argument('--refpdb','--ref', metavar='FILE', help='Reference structure of the protein')
    p.add_argument('--cif', metavar='FILE', help='Other cif files for the structure')
    # Data Arguments
    p.add_argument('--inputmtz','--mtz', metavar='FILE', required=True, help='MTZ reflection data file')
    # Program Arguments
    p.add_argument('--builder', metavar='PROG', default=ligand_builders[0], choices=ligand_builders, help='Program used to generate ligands')
    p.add_argument('--fitter', metavar='PROG', default=ligand_fitters[0], choices=ligand_fitters, help='Program used to fit ligands')
    p.add_argument('--refiner', metavar='PROG', default=refiners[0], choices=refiners, help='Program used to refine structures')
    # File Arguments
    p.add_argument('--outdir', '-o', metavar='DIR', default='./LIGANDMODEL', help='Output Directory')
    p.add_argument('--outtemplate', metavar='DIR', default='', help='Output File Template')
    # Meta
    p.add_argument('--experiment', '--exp', metavar='NAME', default='MODEL', help ='Name describing the experiment i.e. PROTEIN_MODEL_1')
    # Switch Flags
    p.add_argument('--nosubst', action='store_true', default=False, help='Don\'t perform molecular substitution (Dry Run)')
    p.add_argument('--nobuild', action='store_true', default=False, help='Don\'t generate ligand restraints (Dry Run)')
    p.add_argument('--nofit', action='store_true', default=False, help='Don\'t fit ligand (Dry Run)')
    p.add_argument('--norefine', action='store_true', default=False, help='Don\'t refine the fitted protein-ligand complex (Dry Run)')
    p.add_argument('--noscore', action='store_true', default=False, help='Don\'t score the fitted protein-ligand complex (Dry Run)')

    return p


#! /usr/local/python/python2.7.3-64bit/bin/python

# Program Utils
from Bamboo.Utils.Constants import LigandNames
from Bamboo.Utils.ProgramSelector import get_ligand_builder, get_ligand_fitter, get_refiner
# Creating Maps
from Bamboo.Maps.Convert import convert_mtz_to_map
from Bamboo.Maps.Masking import create_masked_map
# Scoring Density
from Bamboo.Density.Edstats.Score import score_against_density_by_subset_file
# Ligand Utils
from Bamboo.Rdkit.Mol import check_ligand_from_file_identical_to_smile, check_ligands_from_file_identical
from Bamboo.Rdkit.Smile import check_smile_readable
from Bamboo.Wrappers.PdbUtils import isolate_residue_by_res_id, remove_residue_by_res_id, merge_pdb_files, reset_pdb_file
from Bamboo.Wrappers.CifUtils import merge_cif_libraries
from Bamboo.Macro.Utils import get_mean_occupancy, get_residue_labels
from Bamboo.Rdkit.Bonds.Fragment import break_and_rename_mol_to_file
# Prototype Pipeline Things
from Prototypes.AutoModel.Utils import modelLog

import os, time

class ligandPipeline(object):
    """Class for automatically modelling Ligands"""

    def __init__(self, inputmtz, refpdb, apopdb, builder, fitter, refiner, outdir, outtemplate='', experiment='ligand-fitting',
                     cif=None, ligsmile=None, ligpdb=None, ligcif=None, verbose=False, **kwargs):

        # Initialise Error Log
        self.error = ErrorObject()
        # Set verboseness
        self.verbose = verbose

        # Check inputs
        if not inputmtz:
            raise SystemExit('No MTZ File supplied to ligand pipeline.')
        if (not builder) and (not kwargs['nobuild']):
            raise SystemExit('No Builder supplied to ligand pipeline.')
        if (not fitter) and (not kwargs['nofit']):
            raise SystemExit('No Fitter supplied to ligand pipeline.')
        if (not refiner) and (not kwargs['nosubst'] or not kwargs['norefine']):
            raise SystemExit('No Refiner supplied to ligand pipeline.')

        print('=======================================================================>')

        # Find what kind of structure the input PDB file is
        if (apopdb and refpdb):
            self.inputpdbtype = 'apo'
            if self.verbose:
                print('==> APO and REF structure supplied:')
                print('=> REF: {!s}'.format(refpdb))
                print('=> APO: {!s}'.format(apopdb))
                print('==> Building from APO structure.')
        elif apopdb:
            self.inputpdbtype = 'apo'
            if self.verbose:
                print('==> APO structure supplied: {!s}'.format(apopdb))
            raise SystemExit('TODO TODO TODO')
        elif refpdb:
            self.inputpdbtype = 'ref'
            if self.verbose:
                print('==> REF structure supplied: {!s}'.format(refpdb))
        else:
            raise SystemExit('Neither APO nor REF structure supplied.')

        # Find what kind of ligand information has been given
        if (ligpdb and ligcif and ligsmile):
            self.ligsmilemol = check_smile_readable(ligsmile)
            self.inputligtype = 'cif'
            if self.verbose:
                print('==> Ligand CIF and SMILE information supplied:')
                print('=> CIF: {!s}'.format(ligcif))
                print('=> SMILE: {!s}'.format(ligsmile))
                print('==> Building using CIF information.')
        elif (ligpdb and ligcif):
            self.inputligtype = 'cif'
            if self.verbose:
                print('==> Ligand CIF information supplied: {!s}'.format(ligcif))
        elif ligsmile:
            self.ligsmilemol = check_smile_readable(ligsmile)
            self.inputligtype = 'smile'
            if self.verbose:
                print('==> Ligand SMILE information supplied: {!s}'.format(ligsmile))
        else:
            raise SystemExit('Neither SMILE nor CIF given for ligand.')

        # Store all inputs (even if None - they will be populated later)
        self.refpdb   = os.path.abspath(refpdb)   if refpdb else None
        self.apopdb   = os.path.abspath(apopdb)   if apopdb else None
        self.cif      = os.path.abspath(cif)      if cif else None
        self.inputmtz = os.path.abspath(inputmtz) if inputmtz else None
        # Ligand Information
        self.ligsmile = ligsmile
        self.ligpdb   = os.path.abspath(ligpdb)   if ligpdb else None
        self.ligcif   = os.path.abspath(ligcif)   if ligcif else None
        # Ligand when fitted to structure (ligand only)
        self.fittedlig = None
        # Fitted ligand added to apo structure
        self.mergedpdb = None
        # Fitter Structure with refined B-factors and Occupancy
        self.factorpdb = None
        self.factormtz = None
        self.factorlog = None
        # Fitted structure after Full Refinement
        self.refinepdb = None
        self.refinemtz = None
        self.refinelig = None
        # Ligand Scores
        self.ligscores = {}
        self.ligscoresfiles = {}
        self.scoresummary = {}

        # Programs to be used
        self.builder = get_ligand_builder(builder)
        self.fitter = get_ligand_fitter(fitter)
        self.refiner = get_refiner(refiner)

        # Set verboseness
        self.builder.verbose = self.verbose
        self.fitter.verbose  = self.verbose
        self.refiner.verbose = self.verbose

        # Output directory
        self.outdir = os.path.abspath(outdir)
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        # Create output directory for the refined APO structure
        self.refiner.apodir = os.path.join(self.outdir,'1-apo')
        if not os.path.exists(self.refiner.apodir):
            os.mkdir(self.refiner.apodir)
        # Create output directory for the builder - LIGAND CIF + PDB
        self.builder.outdir = os.path.join(self.outdir,'2-ligand')
        if not os.path.exists(self.builder.outdir):
            os.mkdir(self.builder.outdir)
        # Create output directory for the fitter - FITTED LIGAND
        self.fitter.outdir = os.path.join(self.outdir,'3-fitted')
        if not os.path.exists(self.fitter.outdir):
            os.mkdir(self.fitter.outdir)
        # Create output directory for the refined structures
        self.refiner.refdir = os.path.join(self.outdir,'4-refined')
        if not os.path.exists(self.refiner.refdir):
            os.mkdir(self.refiner.refdir)

        # Output Scores and Error Log
        self.scorelog = os.path.join(self.outdir,''.join([fitter,'.scores']))
        # TODO TODO TODO Make this log more useful - pass as argument to functions and have all things added to this?
        self.warninglog = os.path.join(self.outdir,''.join(['-'.join(['warnings',builder,fitter]),'.log']))
        # Save the Model log
        self.log = modelLog(ligsmile=self.ligsmile, refpdb=self.refpdb, apopdb=self.apopdb, inputmtz=self.inputmtz, \
                                builder=builder, fitter=fitter, refiner=refiner, outdir=self.outdir, name='ModelLog', verbose=verbose)

        # Create symlink to Reference PDB or APO PDB if present
        if self.refpdb and os.path.exists(self.refpdb) and os.path.exists(self.inputmtz):
            symlink_pdb = os.path.join(self.outdir,'reference.pdb')
            if not os.path.exists(symlink_pdb):
                os.symlink(os.path.relpath(self.refpdb, start=self.outdir), symlink_pdb)
                if self.inputpdbtype == 'ref':
                    os.symlink(os.path.relpath(self.inputmtz, start=self.outdir), os.path.join(self.outdir,'rawdata.mtz'))
        if self.apopdb and os.path.exists(self.apopdb) and os.path.exists(self.inputmtz):
            symlink_pdb = os.path.join(self.outdir, 'apo.pdb')
            if not os.path.exists(symlink_pdb):
                os.symlink(os.path.relpath(self.apopdb, start=self.outdir), symlink_pdb)
                if self.inputpdbtype == 'apo':
                    os.symlink(os.path.relpath(self.inputmtz, start=self.outdir), os.path.join(self.outdir,'apo.mtz'))
        if self.cif and os.path.exists(self.cif):
            symlink_cif = os.path.join(self.outdir,'input_restraints.cif')
            if not os.path.exists(symlink_cif):
                os.symlink(os.path.relpath(self.cif, start=self.outdir), symlink_cif)

        # Output File Template
        self.outtemplate = outtemplate
        # Experiment Label
        self.experiment = experiment

        if self.verbose:
            print('========================================>')
            print('=> {!s:<12} - {!s}'.format('builder', self.builder.name))
            print('=> {!s:<12} - {!s}'.format('fitter', self.fitter.name))
            print('=> {!s:<12} - {!s}'.format('refiner', self.refiner.name))
            print('========================================>')
            print('==> Writing Output Files to: {!s}'.format(self.outdir))
            print('========================================>')

    # ==================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    def PerformMolecularSubstitution(self, refpdb=None, rawmtz=None, cif=None, outdir=None, outfile=None, reset=True, flags=[]):
        """Takes refpdb and refines it against rawmtz - assumes almost the same structure as refpdb was originally solved"""

        # Fill in any fields that have not been provided with the defaults
        if not outdir:
            outdir = self.refiner.apodir
        if not outfile:
            if self.outtemplate:
                outfile = '-'.join(['apo',self.outtemplate,self.refiner.name])
            else:
                outfile = '-'.join(['apo',self.refiner.name])
        if not refpdb:
            refpdb = self.refpdb
        if not rawmtz:
            rawmtz = self.inputmtz
        if not cif:
            cif = self.cif

        if self.verbose:
            print('=> Performing molecular substitution. Flags: {!s}'.format(', '.join([str(f) for f in flags]) if flags else None))

        if reset:
            pdbin = os.path.join(outdir, 'reference.reset.pdb')
            if not os.path.exists(pdbin):
                pdbin = reset_pdb_file(refpdb, pdbin)
        else:
            pdbin = refpdb

        # Pass the command to the Refiner @ self.refiner
        self.apopdb, self.apomtz, self.apolog = self.refiner.run_refinement(inpdb=pdbin, inmtz=rawmtz, incif=cif, outdir=outdir, outfile=outfile, flags=flags)

        # Check for errors
        if hasattr(self.refiner, 'Refiner') and self.refiner.Refiner.err:
            self.error.recordWarning(self.refiner.Refiner.err)
        if hasattr(self.refiner, 'Refiner') and self.refiner.Refiner.timedout:
            raise PipelineError(' {!s} has timed out. See {!s}.'.format(self.refiner.name.upper(), self.apolog))

        if not os.path.exists(self.apopdb):
            raise MolecularSubstitutionError(' {!s} has failed. Output PDB file does not exist. See {!s}'.format(self.refiner.name.upper(), self.apolog))
        if not os.path.exists(self.apomtz):
            raise MolecularSubstitutionError(' {!s} has failed. Output MTZ file does not exist. See {!s}'.format(self.refiner.name.upper(), self.apolog))

        # Everything run smoothly - update variables and add to log
        self.log.setApoPDB(self.apopdb)
        self.log.setApoMTZ(self.apomtz)

        # Create a map of the apo mtz
        Mapper = convert_mtz_to_map(self.apomtz, maptype='FOFC')

    # ==================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    def GenerateLigand(self, ligsmile, outdir=None, outfile=None, flags=[]):
        """Build the Ligand with variable settings - shortcuts to the Builder Method"""

        if outdir == None:
            outdir = self.builder.outdir
        if outfile == None:
            if self.outtemplate:
                outfile = '-'.join(['ligand',self.outtemplate,self.builder.name])
            else:
                outfile = '-'.join(['ligand',self.builder.name])

        if self.verbose:
            print('=> Generating ligand model from SMILE string. Flags: {!s}'.format(', '.join([str(f) for f in flags]) if flags else None))

        # Pass the command to the Builder @ self.builder
        self.ligpdb, self.ligcif, self.liglog = self.builder.generate_ligand(ligsmile=ligsmile, outdir=outdir, outfile=outfile, flags=flags)

        # Check for errors
        if hasattr(self.builder, 'Builder') and self.builder.Builder.err:
            self.error.recordWarning(self.builder.Builder.err)
        if hasattr(self.builder, 'Builder') and self.builder.Builder.timedout:
            raise PipelineError(' {!s} has timed out. See {!s}.'.format(self.builder.name.upper(), self.liglog))

        if not os.path.exists(self.ligpdb):
            raise PipelineError(' {!s} has failed. Output LIG PDB file does not exist. See {!s}'.format(self.builder.name.upper(), self.liglog))
        if not os.path.exists(self.ligcif):
            raise PipelineError(' {!s} has failed. Output LIG CIF file does not exist. See {!s}'.format(self.builder.name.upper(), self.liglog))

        # Check the new file against the smile string
        try:
            rc, message = check_ligand_from_file_identical_to_smile(self.ligsmile, self.ligpdb)
        except LigandCheckError as err:
            raise LigandCheckError(' {!s} has failed. {!s}'.format(self.builder.name.upper(), str(err)))
        if rc:
            raise LigandCheckError(' {!s} has failed. Generated Ligand is invalid: {!s} ({!s})'.format(self.builder.name.upper(), message, self.ligpdb))

        # Add the new cif file to the list of cifs
        if not self.cif:
            self.cif = self.ligcif
        else:
            self.cif = merge_cif_libraries(incifs=[self.cif, self.ligcif], outcif=self.ligcif.replace('.cif','.combined.cif'))

        # Everything run smoothly - update variables and add to log
        self.log.lig.unfitted    = self.ligpdb
        self.log.lig.restraints  = self.ligcif
        self.log.lig.buildinglog = self.liglog

    # ==================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    def FitLigand(self, ligpdb, ligcif, mtz, apopdb, outdir=None, outfile=None, flags=[]):
        """Fit the Ligand with variable settings - shortcuts to the Fitter Method"""

        if outdir == None:
            outdir = self.fitter.outdir
        if outfile == None:
            if self.outtemplate:
                outfile = '-'.join(['fitted',self.outtemplate,self.fitter.name])
            else:
                outfile = '-'.join(['fitted',self.fitter.name])

        if self.verbose:
            print('=> Fitting ligand to APO structure. Flags: {!s}'.format(', '.join([str(f) for f in flags]) if flags else None))

        # Pass the command to the Fitter @ self.fitter
        self.fittedlig, self.mergedpdb, self.fitlog = self.fitter.fit_ligand(ligcif=ligcif, ligpdb=ligpdb, mtz=mtz, apopdb=apopdb, outdir=outdir, outfile=outfile, flags=flags)

        # Get the ligand id (residue number, etc)
        self.ligandid = get_residue_labels(self.fittedlig)[0]

        # Check for errors
        if hasattr(self.fitter, 'Fitter') and self.fitter.Fitter.err:
            self.error.recordWarning(self.fitter.Fitter.err)
        if hasattr(self.fitter, 'Fitter') and self.fitter.Fitter.timedout:
            raise PipelineError(' {!s} has timed out. See {!s}.'.format(self.fitter.name.upper(), self.fitlog))

        if not os.path.exists(self.fittedlig):
            raise PipelineError(' {!s} has failed. Output FITTED LIG file does not exist. See {!s}'.format(self.fitter.name.upper(), self.fitlog))
        if not os.path.exists(self.mergedpdb):
            raise PipelineError(' {!s} has failed. Output MERGED PDB file does not exist. See {!s}'.format(self.fitter.name.upper(), self.fitlog))

        # Check the new file against the template ligand
        try:
            rc, message = check_ligands_from_file_identical(self.ligpdb, self.fittedlig)
        except LigandCheckError as err:
            raise LigandCheckError(' {!s} has failed. {!s}'.format(self.fitter.name.upper(), str(err)))
        if rc:
            raise LigandCheckError(' {!s} has failed. Fitted Ligand is invalid: {!s} ({!s})'.format(self.fitter.name.upper(), message, self.fittedlig))

        # Fragment Ligand
        oldfile = self.fittedlig
        dummy, self.fittedfrags, self.fittedfragscomplex = self.FragmentLigand(self.mergedpdb)
        assert dummy==self.fittedlig

        for maptype in ['FOFC','2FOFC']:
            # Create a masked map of the ligand
            Mapper, Masker = create_masked_map(self.apomtz, self.fittedlig, maptype=maptype)
            # Record Warnings
            if Mapper and Mapper.err:
                self.error.recordWarning(Mapper.err)
            if Masker and Masker.err:
                self.error.recordWarning(Masker.err)

        # Everything run smoothly - update variables and add to log
        self.log.addNewPDB(self.mergedpdb, tag='Fitted01')
        self.log.addNewMTZ(self.log.mtz.current.path, tag='Fitted01')
        self.log.lig.fittinglog = self.fitlog
        self.log.lig.fitted     = self.fittedlig
        self.log.lig.merged     = self.mergedpdb

    # ==================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    def RefineBFactors(self, pdb, mtz, cif=None, outdir=None, outfile=None, flags=[]):
        """Refines ONLY the B-factors of a structure"""

        # Fill in any fields that have not been provided with the defaults
        if not outdir:
            outdir = self.fitter.outdir
        if not outfile:
            if self.outtemplate:
                outfile = '-'.join(['fitted','B',self.outtemplate,self.fitter.name])
            else:
                outfile = '-'.join(['fitted','B',self.fitter.name])
        if not cif:
            cif = self.cif

        # Refine B-factors only
        if 'bonly' not in flags:
            flags.append('bonly')

        if self.verbose:
            print('=> Refining B-factors. Flags: {!s}'.format(', '.join([str(f) for f in flags]) if flags else None))

        # Pass the command to the Refiner @ self.refiner
        self.factorpdb, self.factormtz, self.factorlog = self.refiner.run_refinement(inpdb=pdb, inmtz=mtz, incif=cif, outdir=outdir, outfile=outfile, flags=flags)

        # Check for errors
        if hasattr(self.refiner, 'Refiner') and self.refiner.Refiner.err:
            self.error.recordWarning(self.refiner.Refiner.err)
        if hasattr(self.refiner, 'Refiner') and self.refiner.Refiner.timedout:
            raise PipelineError(' {!s} has timed out. See {!s}.'.format(self.refiner.name.upper(), self.factorlog))

        if not os.path.exists(self.factorpdb):
            raise PipelineError(' {!s} has failed. Output PDB file does not exist. See {!s}'.format(self.refiner.name.upper(), self.factorlog))
        if not os.path.exists(self.factormtz):
            raise PipelineError(' {!s} has failed. Output MTZ file does not exist. See {!s}'.format(self.refiner.name.upper(), self.factorlog))

        # Create a map of the refined mtz
        Mapper = convert_mtz_to_map(self.factormtz, maptype='2FOFC')

        # Fragment Ligand Model
        self.factorlig, self.factorfrags, self.factorfragscomplex = self.FragmentLigand(self.factorpdb)

        # Everything run smoothly - update variables and add to log
        self.log.addNewPDB(self.factorpdb, tag='Refine01-Bfactors')
        self.log.addNewMTZ(self.factormtz, tag='Refine01-Bfactors')
        self.log.lig.factorlog = self.factorlog
        self.log.lig.factor    = self.factorlig

    # ==================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    def RefineComplex(self, pdb, mtz, cif=None, outdir=None, outfile=None, flags=[]):
        """Refine Protein-Ligand Complex"""

        # Fill in any fields that have not been provided with the defaults
        if not outdir:
            outdir = self.refiner.refdir
        if not outfile:
            if self.outtemplate:
                outfile = '-'.join(['refined_complex',self.outtemplate,self.fitter.name])
            else:
                outfile = '-'.join(['refined_complex',self.fitter.name])
        if not cif:
            cif = self.cif

        if self.verbose:
            print('=> Refining Complex. Flags: {!s}'.format(', '.join([str(f) for f in flags]) if flags else None))

        # Pass the command to the Refiner @ self.refiner
        self.refinepdb, self.refinemtz, self.refinelog = self.refiner.run_refinement(inpdb=pdb, inmtz=mtz, incif=cif, outdir=outdir, outfile=outfile, flags=flags)

        # Check for errors
        if hasattr(self.refiner, 'Refiner') and self.refiner.Refiner.err:
            self.error.recordWarning(self.refiner.Refiner.err)
        if hasattr(self.refiner, 'Refiner') and self.refiner.Refiner.timedout:
            raise PipelineError(' {!s} has timed out. See {!s}.'.format(self.refiner.name.upper(), self.refinelog))

        if not os.path.exists(self.refinepdb):
            raise PipelineError(' {!s} has failed. Output PDB file does not exist. See {!s}'.format(self.refiner.name.upper(), self.refinelog))
        if not os.path.exists(self.refinemtz):
            raise PipelineError(' {!s} has failed. Output MTZ file does not exist. See {!s}'.format(self.refiner.name.upper(), self.refinelog))

        # Create a map of the refined mtz
        Mapper = convert_mtz_to_map(self.refinemtz, maptype='2FOFC')

        # Fragment Ligand
        self.refinelig, self.refinefrags, self.refinefragscomplex = self.FragmentLigand(self.refinepdb)

        # Everything run smoothly - update variables and add to log
        self.log.addNewPDB(self.refinepdb, tag='Refine02')
        self.log.addNewMTZ(self.refinemtz, tag='Refine02')
        self.log.lig.refininglog = self.refinelog
        self.log.lig.refined     = self.refinelig

    # ==================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    def FragmentLigand(self, complex):
        "Take a refined model, remove the ligand, fragment it, and then replace it"""

        # Unpack the residue id
        ligname, ligchain, lignum, ligins = self.ligandid

        # Output Ligand File
        ligand         = complex.replace('.pdb','.lig.pdb')
        fraggedligand  = ligand.replace('.lig.pdb','.fragged.lig.pdb')
        fraggedcomplex = complex.replace('.pdb','.fragged.pdb')

        # Pull out the ligand for fragmentation
        if not os.path.exists(ligand):
            isolater = isolate_residue_by_res_id(inpdb=complex, outpdb=ligand, chain=ligchain, resnum=lignum)

        if not os.path.exists(ligand):
            if isolater and isolater.err:
                # TODO WRITE ISOLATER ERRORS TO self.factorlog?!
                self.error.recordWarning(isolater.err)
            raise PipelineError('Failed to isolate ligand model: {!s}'.format(complex))

        # Check the new file against the template ligand
        try:
            rc, message = check_ligands_from_file_identical(self.ligpdb, ligand)
        except LigandCheckError as err:
            raise LigandCheckError(' {!s} has failed. Ligand validation has failed. {!s} ({!s})'.format(self.fitter.name.upper(), str(err), ligand))
        if rc:
            raise LigandCheckError(' {!s} has failed. Extracted Ligand is invalid: {!s} ({!s})'.format(self.fitter.name.upper(), message, ligand))

        # Fragment the Isolated Ligand on rotatable bonds
        if not os.path.exists(fraggedligand):
            dummy = break_and_rename_mol_to_file(ligand, fraggedligand)

            if not os.path.exists(fraggedligand):
                raise PipelineError('Failed to fragment ligand model: {!s}'.format(ligand))

        # Remove the Unfragmented Ligand and replace it with the Fragmented Ligand
        if not os.path.exists(fraggedcomplex):
            # Create a file to hold the structure with the ligand removed
            temp_apo_structure = fraggedcomplex.replace('.pdb','.temp.pdb')
            # Remove the unfragmented ligand
            remover = remove_residue_by_res_id(inpdb=complex, outpdb=temp_apo_structure, chain=ligchain, resnum=lignum)
            # Put fragged ligand in its place
            merger = merge_pdb_files(pdb1=temp_apo_structure, pdb2=fraggedligand, pdbout=fraggedcomplex)
            # Remove the temporary file
            os.remove(temp_apo_structure)
            # TODO WRITE MERGER ERRORS TO self.factorlog?!

        if not os.path.exists(fraggedcomplex):
            raise PipelineError('Failed to merge fragmented ligand model with protein structure: {!s}'.format(fraggedligand))

        return ligand, fraggedligand, fraggedcomplex

    # ==================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    def ScoreLigandWithEdstats(self, lig, merged, mtz):
        """Score the ligand molecule in `merged` identified by `lig`"""

        if self.verbose:
            print('=> Scoring ligand model wth EDSTATS: {!s}'.format(merged))

        # Score with edstats
        labels, self.scoresummary[merged], self.ligscores[merged] = score_against_density_by_subset_file(mtzpath=mtz, pdbpath=merged, subpdb=lig)

        self.ligscoresfiles[merged] = [lig, mtz]

        if self.verbose:
            for lab in labels:
                print('\tResidue scored: {!s}'.format(lab))

    # ==================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    def WriteLigandEdstatsScores(self, outlog):
        """Write the output scores from EDSTATS"""

        if os.path.exists(outlog):
            if self.verbose:
                print('\tLigand score file already exists. DOING NOTHING.')
        else:
            # Ligand scores stored as:
            # ligscores[structure][label] = scores
            # => to extract: ligscores[file][label] -> dict of scores

            # List of scored files
            files = self.ligscores.keys()
            scorefiles =  self.ligscoresfiles
            # Record Scores
            headers = self.ligscores[files[0]].values()[0].keys()

            # Other Scores (META)
            metaheaders = ['template','experiment','occupancy']

            with open(outlog,'w') as outscores:
                outscores.write(','.join(metaheaders+['ResID']+headers+['FULLPDB','LIGPDB','MTZ'])+'\n')
                for scored_file in files:
                    for label in self.ligscores[scored_file].keys():
                        scores = self.ligscores[scored_file][label]
                        values = [str(scores[col]) for col in headers]
                        metavalues = [self.outtemplate, self.experiment, str(get_mean_occupancy(scored_file))]
                        resid = '_'.join(map(str,label))
                        file = os.path.basename(scored_file)
                        outscores.write(','.join(metavalues+[resid]+values+[scored_file]+scorefiles[scored_file])+'\n')

    # ==================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    def ScoreLigandWithMogul(self, lig, merged):
        """Score the ligand molecule in `merged` identified by `lig`"""

        if self.verbose:
            print('=> Scoring ligand model with MOGUL')

    # ==================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>






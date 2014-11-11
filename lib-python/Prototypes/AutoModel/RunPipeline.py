#! /usr/local/python/python2.7.3-64bit/bin/python

import os, sys, glob, random

# Pipeline modules
from Prototypes.AutoModel.Utils import modelLigandParser
from Prototypes.AutoModel.Pipeline import ligandPipeline

def run_ligand_pipeline(args=None):
    """Run the ligand pipeline on custom arguments. Args can either be passed explicitly to the function as a string or, if None, will be parsed from the command line."""

    # Initialise the Command Line Parser
    ArgParser = modelLigandParser()
    # Parse the Command Line Inputs
    Inputs, InputsDict = ArgParser.parseArguments(args)

    # Initialise the main Ligand-Modelling Object (Designed to work directly from ArgParser)
    Modeller = ligandPipeline(**InputsDict)

    if (Inputs.refpdb) and (not Inputs.nosubst) and (not Modeller.apopdb):

        # Don't exit if new ligand is encountered
        flags = ['noexit']

        # Run Molecular Substitution if a reference structure has been provided
        Modeller.PerformMolecularSubstitution(flags=flags)

    if not Inputs.nobuild:

        # Generate the Ligand File with the defaults
        Modeller.GenerateLigand(ligsmile=Modeller.ligsmile)

    if not Inputs.nofit:

        # Fit Ligand to APO File with the defaults
        Modeller.FitLigand(ligpdb=Modeller.builder.outpdb, ligcif=Modeller.builder.outcif, mtz=Modeller.apomtz, apopdb=Modeller.apopdb)

        # Extract Ligand Information
        chainid = Modeller.ligandid[1] if Modeller.ligandid[1].strip() else 'X'
        resnum = Modeller.ligandid[2]
        inscode = Modeller.ligandid[3]

        # Refine the occupancy of the ligand
        flags = [('occupancy',(chainid,resnum,inscode))]
        # Don't exit if new ligand is encountered
        flags.append('noexit')

        # Refine B-factors and Occupancy
        Modeller.RefineBFactors(pdb=Modeller.mergedpdb, mtz=Modeller.apomtz, flags=flags)

        if not Inputs.noscore:
            # Identify Files for scoring
            fitlig = Modeller.fittedlig
            fitpdb = Modeller.mergedpdb
            fitmtz = Modeller.apomtz

            # Check if already done
            if fitpdb not in Modeller.ligscores:
                Modeller.ScoreLigandWithEdstats(lig=fitlig, merged=fitpdb, mtz=fitmtz)

            # Identify Fragged Files for scoring
            fraglig = Modeller.fittedfrags
            fragpdb = Modeller.fittedfragscomplex
            fragmtz = Modeller.apomtz

            # Check if already done
            if fragpdb not in Modeller.ligscores:
                Modeller.ScoreLigandWithEdstats(lig=fraglig, merged=fragpdb, mtz=fragmtz)

            # TODO Should this be scored against the APO density instead?
            # Identify Files for scoring
            lig = Modeller.factorlig
            pdb = Modeller.factorpdb
            mtz = Modeller.factormtz

            # Check if already done
            if pdb not in Modeller.ligscores:
                Modeller.ScoreLigandWithEdstats(lig=lig, merged=pdb, mtz=mtz)

            # Identify Fragged Files for scoring
            fraglig = Modeller.factorfrags
            fragpdb = Modeller.factorfragscomplex
            fragmtz = Modeller.factormtz

            # Check if already done
            if fragpdb not in Modeller.ligscores:
                Modeller.ScoreLigandWithEdstats(lig=fraglig, merged=fragpdb, mtz=fragmtz)

    if not Inputs.norefine:

        # Extract Ligand Information
        chainid = Modeller.ligandid[1] if Modeller.ligandid[1].strip() else 'X'
        resnum = Modeller.ligandid[2]
        inscode = Modeller.ligandid[3]

        # Refine the occupancy of the ligand
        flags = [('occupancy',(chainid,resnum,inscode))]
        # Don't exit if new ligand is encountered
        flags.append('noexit')

        # Refine the structure again (with settings for refining the complex)
        Modeller.RefineComplex(pdb=Modeller.factorpdb, mtz=Modeller.factormtz, flags=flags)

        if not Inputs.noscore:
            # Identify Files
            lig = Modeller.refinelig
            pdb = Modeller.refinepdb
            mtz = Modeller.refinemtz

            # Check if already done
            if pdb not in Modeller.ligscores:
                Modeller.ScoreLigandWithEdstats(lig=lig, merged=pdb, mtz=mtz)

            # Identify Fragged Files for scoring
            fraglig = Modeller.refinefrags
            fragpdb = Modeller.refinefragscomplex
            fragmtz = Modeller.refinemtz

            # Check if already done
            if fragpdb not in Modeller.ligscores:
                Modeller.ScoreLigandWithEdstats(lig=fraglig, merged=fragpdb, mtz=fragmtz)

    # DO SOME STUFF WITH MOGUL HERE
    # TODO TODO TODO TODO TODO TODO

    # Write outputscores
    if Modeller.ligscores:
        Modeller.WriteLigandEdstatsScores(Modeller.scorelog)

    # REPORT RESULTS
    bldt = Modeller.builder.runtime
    fitt = Modeller.fitter.runtime
    reft = Modeller.refiner.runtime

    if Modeller.verbose:
        if (bldt > 0.0) or (fitt > 0.0) or (reft > 0.0):
            print('============================================================================>')
            print('Timings:')
            print('          : Program    : (mins):(secs)')
            print('-------------------------------------------------------------')
            if (bldt > 0.0):
                print('Builder   : {!s:<10} : {:>6n}:{:<02n}'.format(Modeller.builder.name,bldt//60,int(bldt%60)))
            if (fitt > 0.0):
                print('Fitter    : {!s:<10} : {:>6n}:{:<02n}'.format(Modeller.fitter.name,fitt//60,int(fitt%60)))
            if (reft > 0.0):
                print('Refiner   : {!s:<10} : {:>6n}:{:<02n}'.format(Modeller.refiner.name,reft//60,int(reft%60)))

    ScorestoReport = ['Ra','CCSa','ZDa','ZOa']
    fscores = [Modeller.ligscores[Modeller.factorpdb].items()[0][1][a] for a in ScorestoReport]
    rscores = [Modeller.ligscores[Modeller.refinepdb].items()[0][1][a] for a in ScorestoReport]

    if Modeller.verbose:
        print('=======================================================================>')
        print('Results:')
        print('                :   RSR    :   RSCC   :   RSZD   :   RSZO  ')
        print('----------------:----------:----------:----------:-----------')
        print('Fitted Ligand   : {:>8.3f} : {:>8.3f} : {:>8.3f} : {:>8.3f}'.format(*fscores))
        print('Refined Ligand  : {:>8.3f} : {:>8.3f} : {:>8.3f} : {:>8.3f}'.format(*rscores))
        print('=======================================================================>')

    # Write Errors
    if Modeller.error.warnings:
        Modeller.error.writeWarnings(Modeller.warninglog)

    return Modeller


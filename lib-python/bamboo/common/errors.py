
#:::::::::::::::::::::::::::::::::#
# ############################### #
# ###     Error Functions     ### #
# ############################### #
#:::::::::::::::::::::::::::::::::#

def throw_error(message):
    '''Print error message AND EXIT'''
    print "\n==========>"
    print "Terminal error"
    print "==========>"
    print "{!s}".format(message)
    print "==========>"
    print "QUITTING"
    print "==========>\n"
    raise SystemExit('Error Thrown')

def flag_error(message):
    '''Print error message AND DO NOT EXIT'''

    print "\n==========>"
    print "Non-terminal error"
    print "==========>"
    print "{!s}".format(message)
    print "==========>"
    print "CONTINUING"
    print "==========>\n"

class ErrorObject(object):
    """Object to record and report errors"""
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.errorflags = {}
        self.warnings = []

    def __str__(self):
        """Formats error messages for printing"""
        return self.format_warnings()

    def print_warnings(self):
        if self.warnings:
            print(self.format_warnings())
        else:
            print('\n=====>\n=====> No Warnings!\n=====>\n')

    def write_warnings(self, filename):
        with open(filename,'w') as logfile:
            logfile.write(self.format_warnings())
        print('\n=====>\n=====> {!s} Warning(s) written to {!s}\n=====>\n'.format(len(self.warnings),filename))

    def format_warnings(self):
        """Formats the warning messages for printing"""
        warnings = ['\n=====>\n=====> {!s} Warning(s)\n=====>\n'.format(len(self.warnings))]
        warnings.extend(['\n {!s:>3} : {!s}'.format(i+1,warn) for i,warn in enumerate(self.warnings)])
        return ''.join(warnings)

    def record_warning(self, message):
        """Records any NON-TERMINAL error message"""
        self.warnings.append(message)
        if self.verbose:
            flag_error(message)

    def record_flag_warning(self, flag, message):
        """Records any NON-TERMINAL error message (with a flag)"""
        self.record_warning(message)
        if flag in self.errorflags:
            self.errorflags[flag].append(message)
        else:
            self.errorflags[flag] = [message]

    def record_error(self, message):
        """Raises any TERMINAL error message"""
        self.warnings.append(message)
        self.print_warnings()
        throw_error(message)

#:::::::::::::::::::::::::::::::::#
# ############################### #
# ###      Error Classes      ### #
# ############################### #
#:::::::::::::::::::::::::::::::::#

# General Errors

class ExternalProgramError(Exception):
    pass

# Pipline Errors

class PipelineError(Exception):
    pass

class OptionError(Exception):
    pass

class MolecularSubstitutionError(Exception):
    pass

class BuildingError(Exception):
    pass

class FittingError(Exception):
    pass

class RefinementError(Exception):
    pass

# Structure Errors

class MacroMolError(Exception):
    pass

class FragmentationError(Exception):
    pass

class DistanceError(Exception):
    pass

class ScoringError(Exception):
    pass

class AnalysisError(Exception):
    pass

# Ligand/Smile Errors

class LigandCheckError(Exception):
    pass

class SmileCheckError(Exception):
    pass

# RDkit Errors

class RDkitError(Exception):
    pass

class RDkitReadError(Exception):
    pass

class RDkitFragmentationError(Exception):
    pass

class BondAssignmentError(Exception):
    pass

# TODO Change this to MolEqualityError
class EqualityError(Exception):
    pass

# MTZ Errors

class ReflectionException(Exception):
    pass

class LabelError(Exception):
    pass

# Selection Errors

class RefinerSelectionError(Exception):
    pass

class LigandBuilderSelectionError(Exception):
    pass

class LigandFitterSelectionError(Exception):
    pass

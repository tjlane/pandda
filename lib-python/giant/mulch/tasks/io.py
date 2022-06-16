class TaskReturn(object):

    def __init__(self, 
        output, 
        status,
        ):

        self.output = output
        self.status = status

class TaskReturnStatus(object):

    def __init__(self, 
        success,
        errors = None,
        warnings = None,
        ):

        self.success = success
        self.errors = errors
        self.warnings = warnings

class ModelDataInputOutput(object):

    def __init__(self, 
        input_pdb = None, 
        input_mtz = None,
        output_pdb = None, 
        output_mtz = None,
        ):

        self.input_pdb = input_pdb
        self.input_mtz = input_mtz
        self.output_pdb = output_pdb
        self.output_mtz = output_mtz



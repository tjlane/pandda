import giant.logs as lg
logger = lg.getLogger(__name__)

def run_program(prog):

    logger(prog.as_string())

    result = prog.run()

    if result['exitcode'] != 0:

        logger.heading('{} exited with an error'.format(prog.program))

        logger.subheading('Stdout')
        logger(str(result.stdout))
        
        logger.subheading('Stderr')
        logger(str(result.stderr))

        raise Exception(
            '{} exited with an error'.format(
                prog.program
                )
            )

def raise_missing(self, filepath, result):

    logger.subheading('Stdout')
    logger(str(result.stdout))
    
    logger.subheading('Stderr')
    logger(str(result.stderr))

    raise Exception(
        'Failed: {} does not exist'.format(
            str(filepath)
            )
        )

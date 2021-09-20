
from .base_processors import (
    ProcessorNotAvailable,
    ProcessWrapper,
    ProcessorDict,
    Processor,
    ProcessorJoblib,
    )

try: 
    from .luigi_processor import ProcessorLuigiSGE
except: 
    ProcessorLuigiSGE = ProcessorNotAvailable

basic_processor = Processor()

def wrapper_call(func):
    return func()

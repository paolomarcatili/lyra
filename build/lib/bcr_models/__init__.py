
"""
DOC
"""

__version__ = '1.0'

from bcr_models.modelling import IgChain
from bcr_models.modelling import IgComplex
from bcr_models import db

class BCRBaseError(Exception): pass
class BCRParserError(BCRBaseError): pass
class BCRRuntimeError(BCRBaseError): pass
class BCRModellingError(BCRBaseError): pass
class BCRDatabaseError(BCRBaseError): pass

#
# Special Exceptions
# ------------------

class BCRInsertionsInAlignmentError(BCRBaseError): pass

class BCRNoTemplateError(BCRModellingError):

    def __init__(self, msg, ig_chain=None, csid=None):
        super(BCRNoTemplateError, self).__init__(msg)
        self.ig_chain = ig_chain
        self.csid     = csid

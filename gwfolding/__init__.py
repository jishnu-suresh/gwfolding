"""gwfolding: folding of Stochastic Intermediate Data (SID) to one sidereal day.

This package wraps the existing folding pipeline modules. The scientific
routines are unchanged; only the import paths were made package-relative so the
code can be installed with pip.
"""

from .foldSID import foldSID
from .loadSID import loadSID
from .loadSIDParams import loadSIDParams
from .writeFSID import writeFSID
from . import foldUtils

__all__ = [
    "foldSID",
    "loadSID",
    "loadSIDParams",
    "writeFSID",
    "foldUtils",
]

__version__ = "0.1.0"

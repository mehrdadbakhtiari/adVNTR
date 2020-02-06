
from .hmm import *

try:
    from .base import *
except ImportError:
    raise ImportError('Please check the variable, USE_ENHANCED_HMM')

from .hmm import *
import pyximport; pyximport.install()

try:
    from .base import *
except ImportError:
    raise ImportError('Please install hmm package using setup.py')

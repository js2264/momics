"""
momics
~~~~~~

Cloud-native, TileDB-based multi-omics data format.

:author: Jacques Serizay
:license: CC BY-NC 4.0

"""

from .momics import Momics
from .multirangequery import MultiRangeQuery
from .version import __format_version__, __version__
from . import export

__all__ = ["__version__", "__format_version__", "Momics", "MultiRangeQuery", "export"]

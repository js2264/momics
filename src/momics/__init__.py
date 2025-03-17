"""
momics
~~~~~~

Cloud-native, TileDB-based multi-omics data format for machine learning applications.

:author: Jacques Serizay
:license: CC BY-NC 4.0

"""

from . import momics
from . import query
from . import streamer
from . import config
from . import utils
from . import dataset
from .version import __format_version__, __version__

__all__ = [
    "__format_version__",
    "__version__",
    "config",
    "momics",
    "query",
    "streamer",
    "utils",
    "dataset",
]

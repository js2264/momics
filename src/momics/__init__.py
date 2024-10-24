"""
momics
~~~~~~

Cloud-native, TileDB-based multi-omics data format.

:author: Jacques Serizay
:license: CC BY-NC 4.0

"""

from . import momics
from . import query
from . import streamer
from . import config
from . import utils
from .version import __format_version__, __version__

__all__ = [
    "momics",
    "query",
    "streamer",
    "config",
    "utils",
    "__format_version__",
    "__version__",
]

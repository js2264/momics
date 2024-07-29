try:
    from importlib.metadata import PackageNotFoundError, version
except ImportError:
    from importlib_metadata import PackageNotFoundError, version

try:
    __version__ = version("momics")
except PackageNotFoundError:
    __version__ = "unknown"

__format_version__ = 1

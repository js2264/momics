import click

from ..version import __version__

CONTEXT_SETTINGS = {"help_option_names": ["--help"]}


class UnsortedGroup(click.Group):
    """A click Group that lists commands in the order they were added."""

    def list_commands(self, ctx):
        return list(self.commands)


@click.version_option(__version__, "-v", "--version")
@click.group(
    context_settings=CONTEXT_SETTINGS,
    cls=UnsortedGroup,
    epilog="Check out our docs at https://js2264.github.io/momics/ for more details",
)
def cli():
    pass


# Load and register cli subcommands
from . import (
    add,
    create,
    ls,
    query,
)

__all__ = [
    "add",
    "create",
    "ls",
    "query",
]

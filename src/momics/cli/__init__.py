import click

from ..version import __version__

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


class UnsortedGroup(click.Group):
    """A click Group that lists commands in the order they were added."""

    def list_commands(self, ctx):
        return list(self.commands)


@click.version_option(__version__, "-V", "--version")
@click.group(context_settings=CONTEXT_SETTINGS, cls=UnsortedGroup)
def cli():
    """
    Type -h or --help after any subcommand for more information.

    """


# Load and register cli subcommands
from . import (
    dump,
)

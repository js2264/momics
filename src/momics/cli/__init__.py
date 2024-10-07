import click

from ..version import __version__

CONTEXT_SETTINGS = {"help_option_names": ["--help", "-h"]}


class UnsortedGroup(click.Group):
    """A click Group that lists commands in the order they were added."""

    def list_commands(self, ctx):
        return list(self.commands)


@click.version_option(__version__, "-v", "--version")
@click.group(
    context_settings=CONTEXT_SETTINGS,
    cls=UnsortedGroup,
    invoke_without_command=True,
    epilog="Check out our docs at https://js2264.github.io/momics/ for more details",
)
@click.pass_context
def cli(ctx):
    """Command-line software to manage momics repositories."""
    ctx.ensure_object(dict)
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


# Load and register cli subcommands
from . import (
    add,
    binnify,
    cloudconfig,
    create,
    export,
    ls,
    query,
    remove,
)

__all__ = [
    "create",
    "add",
    "ls",
    "query",
    "remove",
    "export",
    "binnify",
    "cloudconfig",
]

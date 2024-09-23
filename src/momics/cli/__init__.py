import click
import multiprocessing

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
    epilog="Check out our docs at https://js2264.github.io/momics/ for more details",
)
@click.option(
    "-@",
    "--threads",
    default=1,
    help="Number of threads to use in parallel operations (default: 1)",
)
@click.pass_context
def cli(ctx, threads):
    """Command-line software to manage momics repositories."""
    ctx.ensure_object(dict)
    ctx.obj["threads"] = threads
    pass


# Load and register cli subcommands
from . import (
    add,
    binnify,
    create,
    export,
    ls,
    query,
    remove,
)

__all__ = [
    "add",
    "create",
    "ls",
    "query",
    "remove",
    "binnify",
    "export",
]

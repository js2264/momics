import click
import cloup
from cloup import Context, Section, Style, HelpFormatter, HelpTheme

from ..version import __version__

CONTEXT_SETTINGS = Context.settings(
    help_option_names=["--help", "-h"],
    # Cloup parameters of Command:
    align_option_groups=False,
    align_sections=True,
    show_constraints=True,
    # Cloup parameters of HelpFormatter:
    formatter_settings=HelpFormatter.settings(
        max_width=110,
        indent_increment=4,
        col_spacing=2,
        theme=HelpTheme.dark().with_(
            heading=Style(fg="white"),
            invoked_command=Style(fg="cyan"),
            col1=Style(fg="cyan"),
        ),
    ),
)


class CloupGroup(cloup.Group):
    """A cloup Group that lists commands in the order they were added."""

    def list_commands(self, ctx):
        return list(self.commands)


class Sections:
    io = Section("Data I/O")
    management = Section("Repository management")
    cloud = Section("Cloud configuration")
    query = Section("Query engine")
    utils = Section("Utils")


@cloup.group(
    context_settings=CONTEXT_SETTINGS,
    cls=CloupGroup,
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
    ingest,
    binnify,
    cloudconfig,
    cp,
    create,
    delete,
    ls,
    query,
    remove,
    tree,
    consolidate,
    manifest,
    version,
)

__all__ = [
    "__version__",
    "binnify",
    "cloudconfig",
    "consolidate",
    "cp",
    "create",
    "delete",
    "ingest",
    "ls",
    "manifest",
    "query",
    "remove",
    "tree",
    "version",
]

if __name__ == "__main__":
    cli()

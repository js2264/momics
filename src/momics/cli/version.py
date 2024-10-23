import click

from ..version import __version__
from .cli import cli
from .cli import Sections


@cli.command(section=Sections.utils)
@click.pass_context
def version(ctx):
    """Print momics version."""
    click.echo(f"momics, version {__version__}")

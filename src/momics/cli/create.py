import click
import cloup

from .. import momics
from .cli import cli
from .cli import Sections


@cli.command(section=Sections.management)
@cloup.argument("path", help="Path to a momics repository", metavar="MOMICS_REPO", required=True)
@click.pass_context
def create(ctx, path):
    """Initiate a Momics repository."""
    path = click.format_filename(path)
    m = momics.Momics(path)
    print(m.path)

import click

from .. import momics
from . import cli


@cli.command()
@click.argument("path", metavar="MOMICS_REPO", required=True)
def create(path):
    """Initiate a Momics repository."""
    path = click.format_filename(path)
    print(path)
    print(momics.Momics(path, create=True))

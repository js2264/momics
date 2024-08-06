import click

from .. import api
from . import cli


@cli.command()
@click.argument("path", metavar="MOMICS_REPO", type=click.Path(exists=False))
def create(path):
    """Initiate a Momics repository."""
    path = click.format_filename(path)
    print(path)
    print(api.Momics(path, create=True))

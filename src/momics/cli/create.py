import click

from .. import momics
from . import cli


@cli.command()
@click.argument("path", metavar="MOMICS_REPO", required=True)
@click.pass_context
def create(ctx, path):
    """Initiate a Momics repository."""
    path = click.format_filename(path)
    m = momics.Momics(path, create=True)
    print(m.path)

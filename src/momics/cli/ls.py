import click

from .. import api
from . import cli


@cli.command()
@click.argument("path", metavar="MOMICS_REPO", type=click.Path(exists=True))
@click.option(
    "--table",
    "-t",
    help="Which supporting table to list.",
    type=click.Choice(["tracks", "chroms"]),
    default="tracks",
    show_default=True,
)
def ls(path, table):
    """List tracks/chromosomes registered in a Momics."""
    if table == "tracks":
        print(api.Momics(path, create=False).tracks())
    if table == "chroms":
        print(api.Momics(path, create=False).chroms())

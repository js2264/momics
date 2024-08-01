import click

from .. import api
from . import cli


@cli.command()
@click.argument("momics", metavar="MOMICS_PATH")
@click.option(
    "--table",
    "-t",
    help="Which table to dump.",
    type=click.Choice(["chroms", "tracks"]),
    default="chroms",
    show_default=True,
)
def dump(momics, table):
    if table == "tracks":
        print(api.Momics(momics).tracks())
    if table == "chroms":
        print(api.Momics(momics).chroms())

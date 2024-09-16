import click
import numpy as np

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
        tr = api.Momics(path, create=False).tracks()
        print(tr.iloc[np.where(tr["label"] != "None")].iloc[:, 0:2])
    if table == "chroms":
        print(api.Momics(path, create=False).chroms())

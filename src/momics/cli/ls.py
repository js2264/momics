import click
import numpy as np

from .. import momics
from . import cli


@cli.command()
@click.argument("path", metavar="MOMICS_REPO", required=True)
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
        tr = momics.Momics(path, create=False).tracks()
        print(
            tr.iloc[np.where(tr["label"] != "None")]
            .iloc[:, 0:2]
            .to_csv(sep="\t", index=False)
        )
    if table == "chroms":
        res = momics.Momics(path, create=False).chroms()
        print(res.iloc[:, 1:].to_csv(sep="\t", index=False, header=False))

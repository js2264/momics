import click
import cloup
import numpy as np

from .. import momics
from .cli import cli
from .cli import Sections


@cli.command(section=Sections.management)
@cloup.argument("path", help="Path to a momics repository", metavar="MOMICS_REPO", required=True)
@click.option(
    "--table",
    "-t",
    help="Which supporting table to list.",
    type=click.Choice(["tracks", "chroms", "features"]),
    default="tracks",
    show_default=True,
)
@click.pass_context
def ls(ctx, path, table):
    """List tracks/chromosomes/features registered in a Momics."""
    if table == "tracks":
        tr = momics.Momics(path).tracks()
        print(tr.iloc[np.where(tr["label"] != "None")].iloc[:, 0:2].to_csv(sep="\t", index=False))
    if table == "features":
        tr = momics.Momics(path).features()
        print(tr.iloc[np.where(tr["label"] != "None")].iloc[:, 0:3].to_csv(sep="\t", index=False))
    if table == "chroms":
        res = momics.Momics(path).chroms()
        print(res.iloc[:, 1:].to_csv(sep="\t", index=False, header=False))

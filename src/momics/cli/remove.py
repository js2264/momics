import click
import numpy as np

from .. import momics
from . import cli


@cli.command()
@click.option(
    "--track",
    "-t",
    help="Track label",
    type=str,
    multiple=True,
    required=True,
)
@click.argument("path", metavar="MOMICS_REPO", required=True)
@click.pass_context
def remove(ctx, path, track):
    """Remove tracks from a momics repo."""
    m = momics.Momics(path)
    for tr in track:
        m.remove_track(tr)
    print(m.tracks().iloc[np.where(m.tracks()["label"] != "None")].iloc[:, 0:2])

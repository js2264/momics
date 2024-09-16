import click
import numpy as np

from .. import api
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
@click.argument("path", metavar="MOMICS_REPO", type=click.Path(exists=True))
def remove(path, track):
    """Remove tracks from a momics repo."""
    m = api.Momics(path, create=False)
    for tr in track:
        m.remove_track(tr)
    print(m.tracks().iloc[np.where(m.tracks()["label"] != "None")].iloc[:, 0:2])

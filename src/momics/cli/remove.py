import click
import cloup
import numpy as np

from .. import momics
from .cli import cli
from .cli import Sections


@cli.command(section=Sections.io)
@click.option(
    "--track",
    "-t",
    help="Track label",
    type=str,
    multiple=True,
    required=True,
)
@cloup.argument("path", help="Path to a momics repository", metavar="MOMICS_REPO", required=True)
@click.pass_context
def remove(ctx, path, track):
    """Remove tracks from a momics repo."""
    m = momics.Momics(path)
    for tr in track:
        m.remove_track(tr)
    print(m.tracks().iloc[np.where(m.tracks()["label"] != "None")].iloc[:, 0:2])

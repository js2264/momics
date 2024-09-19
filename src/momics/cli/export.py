import click
import numpy as np

from .. import api
from . import cli


@cli.command()
@click.option(
    "--track",
    "-t",
    help="Track to export as bigwig",
    type=str,
    required=True,
)
@click.option(
    "--prefix",
    "-p",
    help="Prefix of bigwig file to write",
    type=str,
    required=True,
)
@click.argument("path", metavar="MOMICS_REPO", type=click.Path(exists=True))
def export(path, track, prefix):
    """Export a track from a momics repo as a bigwig file."""
    m = api.Momics(path, create=False)
    m.export_track(track, f"{prefix}.bw")

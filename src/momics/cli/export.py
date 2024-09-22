import click

from .. import momics
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
    "--output",
    "-o",
    help="Path of bigwig file to write",
    type=str,
    required=True,
)
@click.argument("path", metavar="MOMICS_REPO", required=True)
def export(path, track, output):
    """Export a track from a momics repo as a bigwig file."""
    m = momics.Momics(path, create=False)
    m.export_track(track, output)

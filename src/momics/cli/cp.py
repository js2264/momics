import click

from .. import export, momics
from . import cli


@cli.command()
@click.option(
    "--track",
    "-t",
    help="Track to copy as a bigwig file",
    type=str,
    required=False,
)
@click.option(
    "--features",
    "-f",
    help="Features set to copy as a bed file",
    type=str,
    required=False,
)
@click.option(
    "--output",
    "-o",
    help="Path of output file to write",
    type=str,
    required=True,
)
@click.argument("path", metavar="MOMICS_REPO", required=True)
@click.pass_context
def cp(ctx, path, track, features, output):
    """Copy a track/feature set from a momics repo to a bigwig/bed file."""
    m = momics.Momics(path)
    if track:
        export.export_track(m, track, output)
    elif features:
        export.export_features(m, features, output)
    else:
        print("Please specify either --track or --features")
        return False

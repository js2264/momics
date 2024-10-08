import os
import click

from .. import export, momics
from . import cli


@cli.command()
@click.option(
    "--type",
    "-t",
    help="Type of data to extract",
    type=click.Choice(["sequence", "track", "features"]),
    required=True,
)
@click.option(
    "--label",
    "-l",
    help="For `track` and `features` types, the name of the track or feature set to extract.",
    type=str,
    required=False,
)
@click.option(
    "--output",
    "-o",
    help="Path of output file to write.",
    type=str,
    required=True,
)
@click.option(
    "--force",
    "-f",
    help="Force overwrite of existing files.",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.argument("path", metavar="MOMICS_REPO", required=True)
@click.pass_context
def cp(ctx, path, type, label, force, output):
    """Copy sequence/track/feature set from a momics repo to a fa/bigwig/bed file."""

    if not force:
        click.confirm(
            f"{output} file already exists. \
            Are you sure you want to overwrite it",
            abort=True,
        )
        os.remove(output)

    m = momics.Momics(path)
    if type == "sequence":
        export.export_sequence(m, output)
    elif type == "track":
        export.export_track(m, label, output)
    elif type == "features":
        export.export_features(m, label, output)
    else:
        return False

    return True

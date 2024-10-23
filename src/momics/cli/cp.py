import os
import click
import cloup

from .. import momics
from .cli import cli
from .cli import Sections


@cli.command(section=Sections.io)
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
@cloup.argument("path", help="Path to a momics repository", metavar="MOMICS_REPO", required=True)
@click.pass_context
def cp(ctx, path, type, label, force, output):
    """Copy sequence/track/feature set from a momics repo to a fa/bigwig/bed file."""

    if not force and output:
        if os.path.exists(output):
            click.confirm(
                f"{output} file already exists. \
                Are you sure you want to overwrite it",
                abort=True,
            )
            os.remove(output)

    m = momics.Momics(path)
    if type == "sequence":
        m.export_sequence(output)
    elif type == "track":
        m.export_track(label, output)
    elif type == "features":
        m.export_features(label, output)
    else:
        return False

    return True

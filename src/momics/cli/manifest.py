import json
import os
import click
import cloup

from momics import momics as m

from .cli import cli
from .cli import Sections


@cli.command(section=Sections.management)
@click.option(
    "--output",
    "-o",
    help="Path of output JSON file to write.",
    type=str,
    required=False,
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
def manifest(ctx, path, output, force):
    """Print the manifest of a momics repository."""

    if not force and output:
        if os.path.exists(output):
            click.confirm(
                f"{output} file already exists. Are you sure you want to overwrite it",
                abort=True,
            )
            os.remove(output)

    mom = m.Momics(path)
    man = mom.manifest()

    if not output:
        print(json.dumps(man, indent=2))
        return None
    else:
        with open(output, "w") as f:
            json.dump(man, f, indent=2)
        return True

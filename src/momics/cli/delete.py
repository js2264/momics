import os

import click
import cloup

from .. import momics
from .cli import cli
from .cli import Sections


@cli.command(section=Sections.management)
@click.option(
    "--yes",
    "-y",
    help="Delete the repository without asking for confirmation.",
    is_flag=True,
    default=False,
    show_default=True,
)
@cloup.argument("path", help="Path to a momics repository", metavar="MOMICS_REPO", required=True)
@click.pass_context
def delete(ctx, path, yes):
    """Delete a momics repo."""
    if not os.path.exists(path):
        print(f"Repository {path} does not exist.")
        return False
    m = momics.Momics(path)
    if not yes:
        click.confirm(f"Are you sure you want to delete the momics repo at {path}?", abort=True)
        m.remove()
    else:
        m.remove()
    return True

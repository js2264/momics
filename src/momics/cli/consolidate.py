import logging
import click
import cloup

from momics import momics as m

from .cli import cli
from .cli import Sections


@cli.command(section=Sections.management)
@click.option(
    "--vacuum",
    "-v",
    help="Flag to indicate wether fragments should be vacuumed after consolidation.",
    is_flag=True,
    required=False,
)
@cloup.argument("path", help="Path to a momics repository", metavar="MOMICS_REPO", required=True)
@click.pass_context
def consolidate(ctx, path, vacuum):
    """Consolidate a momics repository.

    Consolidation is the process of compacting all the fragments from the repository
    to remove any unused space. This process is useful to reduce the size of the
    repository and improve performance.

    The vacuum flag indicates wether the consolidated
    array should be vacuumed after consolidation. This process is useful to further
    reduce the size of the repository.
    """

    mom = m.Momics(path)
    s = mom.size()
    mom.consolidate(vacuum=vacuum)
    sf = mom.size()
    logging.info(f"Final repository size: {sf/(1024*1024)} Mb.")
    logging.info(f"Space saved after consolidation: {(s - sf)/(1024*1024)} Mb.")
    return True

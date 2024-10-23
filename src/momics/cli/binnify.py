import cloup
import click

from .. import momics
from .cli import cli
from .cli import Sections


@cli.command(section=Sections.utils)
@click.option(
    "--width",
    "-w",
    help="The width of each bin.",
    type=int,
    default=None,
    required=True,
)
@click.option(
    "--step",
    "-s",
    help="The step size for tiling.",
    type=int,
    default=None,
    required=True,
)
@click.option(
    "--cut_last_bin_out",
    "-c",
    help="Remove the last bin in each chromosome, which likely does not have " + "the same width.",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    "--output",
    "-o",
    help="Path to a bed file to write",
    type=str,
    required=False,
)
@cloup.argument("path", help="Path to a momics repository", metavar="MOMICS_REPO", required=True)
@click.pass_context
def binnify(ctx, path, width, step, cut_last_bin_out, output):
    """Binnify chromosomes from a Momics repository."""
    m = momics.Momics(path)
    bins = m.bins(width, step, cut_last_bin_out).df
    bins = bins.to_csv(sep="\t", index=False, header=False)

    if output is not None:
        with open(output, "w") as file:
            file.write(bins)
    else:
        print(bins)

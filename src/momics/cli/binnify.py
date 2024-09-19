import click

from .. import api
from . import cli


@cli.command()
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
    help="Remove the last bin in each chromosome, which likely does not have the same width.",
    type=bool,
    default=True,
    show_default=True,
    required=False,
)
@click.option(
    "--prefix",
    "-p",
    help="Prefix of bed file to write",
    type=str,
    required=False,
)
@click.argument("path", metavar="MOMICS_REPO", type=click.Path(exists=False))
def binnify(path, width, step, cut_last_bin_out, prefix):
    """Binnify chromosomes from a Momics repository."""
    m = api.Momics(path, create=False)
    bins = m.bins(width, step, cut_last_bin_out)
    bins = bins.to_csv(sep="\t", index=False, header=False)
    
    if prefix is not None:
        output = f"{prefix}.bed"
        with open(output, "w") as file:
            file.write(bins)
    else:
        print(bins)

import click

from .. import api
from . import cli


@cli.group()
@click.pass_context
def add(ctx):
    """Add a file to a Momics."""


@add.command()
@click.option(
    "--genome",
    "-g",
    help="Genome reference (e.g. hg38, sacCer3, ...).",
    default="",
)
@click.argument(
    "file", metavar="CHROM_SIZES", required=True, type=click.Path(exists=True)
)
@click.argument("path", metavar="MOMICS_REPO", required=True)
@click.pass_context
def chroms(ctx, file, genome, path):
    """Register chromosomes sizes to Momics."""
    m = api.Momics(path, create=False)
    chrom_lengths = {}
    with open(file) as chroms:
        for line in chroms:
            chrom, length = line.strip().split()
            chrom_lengths[chrom] = int(length)
    m.add_chroms(chrom_lengths, genome_version=genome)
    print(m.chroms())


@add.command()
@click.argument("files", metavar="FILE ...", nargs=-1, required=True)
@click.argument(
    "path", metavar="MOMICS_REPO", nargs=1, required=True, type=click.Path(exists=True)
)
@click.pass_context
def tracks(ctx, files, path):
    """Add tracks to Momics."""
    fs = {}
    for f in files:
        fs[f.split("=", 1)[0]] = f.split("=", 1)[1]
    m = api.Momics(path, create=False)
    m.add_tracks(fs)
    print(m.tracks())

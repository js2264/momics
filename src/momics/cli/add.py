import click
import numpy as np

from .. import api
from . import cli


@cli.group()
@click.pass_context
def add(ctx):
    """Add a file to a Momics."""


@add.command()
@click.option(
    "--file",
    "-f",
    help="UCSC-style coordinates",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--genome",
    "-g",
    help="Genome reference (e.g. hg38, sacCer3, ...).",
    default="",
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
@click.option(
    "--file",
    "-f",
    help="Named track file, provided as `--file key=value` "
    + "(e.g. `--file bw1=my_file.bw`). The `--file` option can be provided "
    + "several times.",
    type=str,
    multiple=True,
    required=True,
)
@click.argument("path", metavar="MOMICS_REPO", required=True)
@click.pass_context
def tracks(ctx, file, path):
    """Add tracks to Momics."""
    fs = {}
    for f in file:
        fs[f.split("=", 1)[0]] = f.split("=", 1)[1]
    m = api.Momics(path, create=False)
    m.add_tracks(fs)
    print(m.tracks().iloc[np.where(m.tracks()["label"] != "None")].iloc[:, 0:2])


@add.command()
@click.option(
    "--file",
    "-f",
    help="Fasta file",
    type=click.Path(exists=True),
    required=True,
)
@click.argument("path", metavar="MOMICS_REPO", required=True)
@click.pass_context
def seq(ctx, file, path):
    """Add genomic sequence to Momics."""
    m = api.Momics(path, create=False)
    m.add_sequence(file)
    print(m.seq())

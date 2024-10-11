import click
import numpy as np
import pyranges as pr

from .. import momics
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
    m = momics.Momics(path)
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
@click.option(
    "-@",
    "--threads",
    default=1,
    help="Number of threads to use in parallel operations (default: 1)",
)
@click.pass_context
def tracks(ctx, file, path, threads):
    """Add tracks to Momics."""
    fs = {}
    for f in file:
        fs[f.split("=", 1)[0]] = f.split("=", 1)[1]
    m = momics.Momics(path)
    m.add_tracks(fs, threads=threads)
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
@click.option(
    "-@",
    "--threads",
    default=1,
    help="Number of threads to use in parallel operations (default: 1)",
)
@click.pass_context
def seq(ctx, file, path, threads):
    """Add genomic sequence to Momics."""
    m = momics.Momics(path)
    m.add_sequence(file, threads=threads)
    print(m.seq())


@add.command()
@click.option(
    "--file",
    "-f",
    help="Named BED file, provided as `--file key=value` \
        (e.g. `--file bw1=my_file.bw`). The `--file` option can be provided \
        several times. \
        The first three columns of the BED file must describe the genomic \
        coordinates of the features (chromosome, start, end).",
    type=str,
    multiple=True,
    required=True,
)
@click.argument("path", metavar="MOMICS_REPO", required=True)
@click.option(
    "-@",
    "--threads",
    default=1,
    help="Number of threads to use in parallel operations (default: 1)",
)
@click.pass_context
def features(ctx, file, path, threads):
    """Add genomic features to Momics."""
    fs = {}
    for f in file:
        bed = f.split("=", 1)[1]
        bed = pr.read_bed(bed)
        fs[f.split("=", 1)[0]] = bed
    m = momics.Momics(path)
    m.add_features(fs, threads=threads)
    print(m.features().iloc[np.where(m.features()["label"] != "None")].iloc[:, 0:2])

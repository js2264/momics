import glob
import os
import click
import cloup
import pyBigWig
import numpy as np
import pyranges as pr

from ..momics import Momics
from .cli import cli
from .cli import Sections


@cli.group(section=Sections.io)
@click.pass_context
def ingest(ctx):
    """Ingest data file(s) to a Momics repository."""


@ingest.command()
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
@cloup.argument("path", help="Path to a momics repository", metavar="MOMICS_REPO", required=True)
@click.pass_context
def chroms(ctx, file, genome, path):
    """Register chromosomes sizes to Momics."""
    m = Momics(path)
    chrom_lengths = {}
    with open(file) as chroms:
        for line in chroms:
            chrom, length = line.strip().split()
            chrom_lengths[chrom] = int(length)
    m.ingest_chroms(chrom_lengths, genome_version=genome)
    print(m.chroms())


@ingest.command()
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
@cloup.argument("path", help="Path to a momics repository", metavar="MOMICS_REPO", required=True)
@click.option(
    "-@",
    "--threads",
    default=1,
    help="Number of threads to use in parallel operations (default: 1)",
)
@click.pass_context
def tracks(ctx, file, path, threads):
    """Ingest tracks to Momics."""
    fs = {}
    for f in file:
        fs[f.split("=", 1)[0]] = f.split("=", 1)[1]
    m = Momics(path)
    m.ingest_tracks(fs, threads=threads)
    print(m.tracks().iloc[np.where(m.tracks()["label"] != "None")].iloc[:, 0:2])


@ingest.command()
@click.option(
    "--file",
    "-f",
    help="Fasta file",
    type=click.Path(exists=True),
    required=True,
)
@cloup.argument("path", help="Path to a momics repository", metavar="MOMICS_REPO", required=True)
@click.option(
    "-@",
    "--threads",
    default=1,
    help="Number of threads to use in parallel operations (default: 1)",
)
@click.pass_context
def seq(ctx, file, path, threads):
    """Ingest genomic sequence to Momics."""
    m = Momics(path)
    m.ingest_sequence(file, threads=threads)
    print(m.seq())


@ingest.command()
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
@cloup.argument("path", help="Path to a momics repository", metavar="MOMICS_REPO", required=True)
@click.option(
    "-@",
    "--threads",
    default=1,
    help="Number of threads to use in parallel operations (default: 1)",
)
@click.pass_context
def features(ctx, file, path, threads):
    """Ingest genomic features to Momics."""
    fs = {}
    for f in file:
        bed = f.split("=", 1)[1]
        bed = pr.read_bed(bed)
        fs[f.split("=", 1)[0]] = bed
    m = Momics(path)
    m.ingest_features(fs, threads=threads)
    print(m.features().iloc[np.where(m.features()["label"] != "None")].iloc[:, 0:2])


@ingest.command()
@click.option(
    "--folder",
    "-F",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--genome",
    "-g",
    help="Genome reference (e.g. hg38, sacCer3, ...).",
    default="",
)
@click.option(
    "-@",
    "--threads",
    default=1,
    help="Number of threads to use in parallel operations (default: 1)",
)
@cloup.argument("path", help="Path to a momics repository", metavar="MOMICS_REPO", required=True)
@click.pass_context
def bulk(ctx, folder, genome, threads, path):
    """Bulk ingest all bigwig files from a folder into a Momics repository."""
    m = Momics(path)

    bw_files = glob.glob(os.path.join(folder, "*.bw"))
    fa_file = glob.glob(os.path.join(folder, "*.fa"))
    bed_files = glob.glob(os.path.join(folder, "*.bed"))

    # If more than 1 fasta file found and sequence is not already registered, abort
    needs_fa = True
    if len(fa_file) > 0 and not needs_fa:
        ctx.warn("Fasta file found but no sequence already registered. Ignoring fasta.")
    if len(fa_file) > 1:
        ctx.fail("At most one fasta file can be provided. Aborting now.")

    # Ingest chroms
    if m.chroms().empty:
        with pyBigWig.open(bw_files[0]) as bw:
            chroms = bw.chroms()
        m.ingest_chroms(chroms, genome_version=genome)

    # Ingest seq
    if needs_fa and len(fa_file) > 0:
        m.ingest_sequence(fa_file[0], threads=threads)

    # Ingest tracks
    if len(bw_files) > 0:
        bws = {os.path.basename(f).split(".bw")[0]: f for f in bw_files}
        m.ingest_tracks(bws, threads=threads)

    # Ingest features
    if len(bed_files) > 0:
        beds = {os.path.basename(f).split(".bed")[0]: f for f in bed_files}
        gr = {k: pr.read_bed(v) for k, v in beds.items()}
        m.ingest_features(gr, threads=threads)

import click
from Bio import SeqIO
import pyranges as pr

from momics.momicsquery import MomicsQuery

from .. import momics
from .. import utils
from ..logging import logger
from . import cli


def _validate_exclusive_options(file, coordinates):
    if file and coordinates:
        raise click.BadParameter("You must provide either --file or --coordinates, not both.")
    if not file and not coordinates:
        raise click.BadParameter("You must provide one of --file or --coordinates.")


@cli.group()
@click.pass_context
def query(ctx):
    """Query a Momics table"""
    pass


@query.command()
@click.option(
    "--coordinates",
    "-c",
    help="UCSC-style coordinates",
    type=str,
)
@click.option(
    "--file",
    "-f",
    help="BED file listing coordinates to query. If provided, `coordinates` " + "is ignored.",
    type=click.Path(exists=True),
)
@click.option(
    "--output",
    "-o",
    help="Output file to save data in (data will be exported as tsv)",
    type=click.Path(),
    required=False,
    default=None,
    show_default=True,
)
@click.option(
    "-@",
    "--threads",
    default=1,
    help="Number of threads to use in parallel operations (default: 1)",
)
@click.argument("path", metavar="MOMICS_REPO", required=True)
@click.pass_context
def tracks(ctx, path, coordinates, file, output: str, threads: int = 1):
    """Extract track coverages over a chromosome interval."""

    # Validate that either `file` or `coordinates` is provided, but not both
    _validate_exclusive_options(file, coordinates)

    mom = momics.Momics(path)

    if coordinates is not None:
        bed = utils.parse_ucsc_coordinates(coordinates)
    else:
        bed = pr.read_bed(file)

    res = MomicsQuery(mom, bed).query_tracks(threads=threads).to_df()
    if output is None:
        print(res.to_csv(sep="\t", index=False))
    else:
        logger.info(f"Writing coverage data to {output} file...")
        res.to_csv(path_or_buf=output, sep="\t", index=False)


@query.command()
@click.option(
    "--coordinates",
    "-c",
    help="UCSC-style coordinates",
    type=str,
)
@click.option(
    "--file",
    "-f",
    help="BED file listing coordinates to query. If provided, `coordinates` " + "is ignored.",
    type=click.Path(exists=True),
)
@click.option(
    "--output",
    "-o",
    help="Output file to save data in (data will be exported as fasta)",
    type=click.Path(),
    required=False,
    default=None,
    show_default=True,
)
@click.option(
    "-@",
    "--threads",
    default=1,
    help="Number of threads to use in parallel operations (default: 1)",
)
@click.argument("path", metavar="MOMICS_REPO", required=True)
@click.pass_context
def seq(ctx, path, coordinates, file, output: str, threads: int = 1):
    """Extract chromosomal sequences over chromosome intervals."""

    # Validate that either `file` or `coordinates` is provided, but not both
    _validate_exclusive_options(file, coordinates)

    mom = momics.Momics(path)

    if coordinates is not None:
        bed = utils.parse_ucsc_coordinates(coordinates)
    else:
        bed = pr.read_bed(file)

    res = MomicsQuery(mom, bed).query_sequence(threads=threads).to_SeqRecord()
    if output is None:
        for record in res:
            print(f">{record.id}")
            print(record.seq)
    else:
        logger.info(f"Writing sequences to {output} file...")
        with open(output, "w") as fasta_file:
            SeqIO.write(res, fasta_file, "fasta")

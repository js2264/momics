import click
import typing

from .. import api
from . import cli


@cli.group()
@click.pass_context
def query(ctx):
    """Query a Momics table"""


@query.command()
@click.option(
    "--coordinates",
    "-c",
    help="UCSC-style coordinates",
    type=str,
    required=True,
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
@click.argument("path", metavar="MOMICS_REPO", type=click.Path(exists=True))
@click.pass_context
def tracks(ctx, path, coordinates, output: str):
    """Extract track coverages over a chromosome interval."""
    res = api.Momics(path, create=False).query_tracks(coordinates)
    if output is not None:
        res.to_csv(path_or_buf=output, sep="\t", index=False)
    else:
        print(res)


@query.command()
@click.option(
    "--coordinates",
    "-c",
    help="UCSC-style coordinates",
    type=str,
    required=True,
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
@click.argument("path", metavar="MOMICS_REPO", type=click.Path(exists=True))
@click.pass_context
def seq(ctx, path, coordinates, output: str):
    """Extract chromosomal sequence over a chromosome interval."""
    seq = api.Momics(path, create=False).query_sequence(coordinates)
    seq = "".join(seq)
    if output is not None:
        with open(output, "w") as file:
            file.write(f">{coordinates}\n")
            # Split the sequence into lines of 60 characters
            for i in range(0, len(seq), 60):
                file.write(seq[i : i + 60] + "\n")
    else:
        print(seq)

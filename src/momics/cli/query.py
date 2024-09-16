import click
import typing

from .. import api
from . import cli


@cli.command()
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
def query(path, coordinates, output: str):
    """Extract coverages over a chromosome interval."""
    res = api.Momics(path, create=False).query(coordinates)
    if output is not None:
        res.to_csv(path_or_buf=output, sep="\t", index=False)
    else:
        print(res)

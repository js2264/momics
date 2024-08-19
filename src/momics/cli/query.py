import click

from .. import api
from . import cli


@cli.command()
@click.argument("query", metavar="QUERY", type=str)
@click.argument("path", metavar="MOMICS_REPO", type=click.Path(exists=True))
def query(path, query):
    """Extract coverages over a chromosome interval.

    QUERY should follow UCSC-like genomic coordinate convention, e.g. "II:991-1010"
    """
    print(api.Momics(path, create=False).query(query))

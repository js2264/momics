import click
import typing
import pandas as pd

from momics.multirangequery import MultiRangeQuery

from .. import api
from .. import utils
from . import cli


@cli.group()
@click.pass_context
def query(ctx):
    """Query a Momics table"""


def _validate_exclusive_options(file, coordinates):
    if file and coordinates:
        raise click.BadParameter(
            "You must provide either --file or --coordinates, not both."
        )
    if not file and not coordinates:
        raise click.BadParameter("You must provide one of --file or --coordinates.")


def _parse_dict_to_df(res):
    # Prepare empty long DataFrame without scores, to merge with results
    ranges = list(res[list(res.keys())[0]].keys())
    ranges_str = []
    for i, coords in enumerate(ranges):
        chrom, range_part = coords.split(":")
        start = int(range_part.split("-")[0])
        end = int(range_part.split("-")[1])
        w = end - start + 1
        label = [{"chr": chrom, "position": x} for x in range(start, end + 1)]
        ranges_str.extend(label)
    df = pd.DataFrame(ranges_str)
    for track in list(res.keys()):
        df[track] = [value for sublist in res[track].values() for value in sublist]
    return df


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
    help="BED file listing coordinates to query. If provided, `coordinates` is ignored.",
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
@click.argument("path", metavar="MOMICS_REPO", type=click.Path(exists=True))
@click.pass_context
def tracks(ctx, path, coordinates, file, output: str):
    """Extract track coverages over a chromosome interval."""

    # Validate that either `file` or `coordinates` is provided, but not both
    _validate_exclusive_options(file, coordinates)

    mom = api.Momics(path, create=False)

    if coordinates is not None:
        chr, range_part = coordinates.split(":")
        start = int(range_part.split("-")[0])
        end = int(range_part.split("-")[1])
        bed = pd.DataFrame([{"chr": chr, "start": start, "end": end}])
    else:
        bed = utils.import_bed_file(file)

    res = MultiRangeQuery(mom, bed).query_tracks().to_df()
    if output is None:
        print(res.to_csv(sep="\t", index=False))
    else:
        print(res)
        res.to_csv(path_or_buf=output, sep="\t", index=False)


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

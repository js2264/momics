import click
import cloup
import tiledb
import json

from .. import momics
from .cli import cli
from .cli import Sections


@cli.command(section=Sections.utils)
@cloup.argument("path", help="Path to a momics repository", metavar="MOMICS_REPO", required=True)
@click.pass_context
def info(ctx, path):
    """Get info about tracks, features, or chromosomes in a Momics repository."""
    mom = momics.Momics(path)
    chrs = mom.chroms()["chrom"]

    if chrs.empty:
        has_seq = False
        timestamp = "unknown"
        genome_assembly = "unknown"
        genome_length = 0
        has_seq = False
        nchroms = 0
        ntracks = 0
        nfeatures = 0

    else:
        tdb = mom._build_uri("genome", "chroms.tdb")
        with tiledb.open(tdb, "r", ctx=mom.cfg.ctx) as A:
            timestamp = A.meta["timestamp"]
            genome_assembly = A.meta["genome_assembly_version"]
            genome_length = sum(A.df[:]["length"])
            nchroms = A.df[:].shape[0]

        vfs = mom.cfg.vfs
        sequence_uri = mom._build_uri("genome", chrs[0]) + ".tdb"
        tracks_uri = mom._build_uri("coverage", "tracks") + ".tdb"
        features_uri = mom._build_uri("annotations", "features") + ".tdb"
        has_seq = vfs.is_dir(sequence_uri)
        if vfs.is_dir(tracks_uri):
            ntracks = mom.tracks().shape[0]
        else:
            ntracks = 0

        if vfs.is_dir(features_uri):
            nfeatures = mom.features().shape[0]
        else:
            nfeatures = 0

    info = {
        "creation-date": timestamp,
        "format": "TileDB::Momics",
        "format-url": "https://github.com/js2264/momics",
        "genome-assembly": genome_assembly,
        "genome-length": genome_length,
        "genome-sequence": has_seq,
        "nchroms": nchroms,
        "ntracks": ntracks,
        "nfeatures": nfeatures,
    }
    click.echo(json.dumps(info, sort_keys=False, indent=4, separators=(",", ": ")))

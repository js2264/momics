import click

import momics

from . import cli


@cli.command()
@click.argument("path", metavar="MOMICS_REPO", required=True)
@click.pass_context
def tree(ctx, path):
    """List all the TileDB tables already ingested."""

    mom = momics.Momics(path)
    vfs = mom.cfg.vfs
    genome_uri = mom._build_uri("genome")
    chroms_uri = mom._build_uri("genome", "chroms") + ".tdb"
    sequence_uri = mom._build_uri("genome", "sequence")
    tracks_uri = mom._build_uri("coverage", "tracks") + ".tdb"
    features_uri = mom._build_uri("features", "features") + ".tdb"
    print(f"Momics repository    {path}")
    print("-----------------    " + "".join(["-"] * len(path)))
    if vfs.is_dir(chroms_uri):
        print(f"|- genome            {genome_uri}")
        print(f"   |- chroms         {chroms_uri}")
    if vfs.is_dir(sequence_uri):
        print(f"   |- sequence       {sequence_uri}")
    if vfs.is_dir(tracks_uri):
        print(f"|- tracks            {tracks_uri}")
    if vfs.is_dir(features_uri):
        print(f"|- features          {features_uri}")
    return True

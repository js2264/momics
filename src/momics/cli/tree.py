import os
from pathlib import Path
import click

import momics

from . import cli


@cli.command()
@click.argument("path", metavar="MOMICS_REPO", required=True)
@click.pass_context
def tree(ctx, path):
    """List all the TileDB tables already ingested."""

    mom = momics.Momics(path)
    chrs = mom.chroms()["chrom"]
    name = Path(os.path.dirname(path)).with_suffix("").name
    vfs = mom.cfg.vfs
    chroms_uri = mom._build_uri("genome", "chroms") + ".tdb"
    sequence_uri = mom._build_uri("genome", chrs[0]) + ".tdb"
    tracks_uri = mom._build_uri("coverage", "tracks") + ".tdb"
    features_uri = mom._build_uri("annotations", "features") + ".tdb"

    has_chroms = vfs.is_dir(chroms_uri)
    has_seq = vfs.is_dir(sequence_uri)
    has_tracks = vfs.is_dir(tracks_uri)
    has_features = vfs.is_dir(features_uri)

    if has_chroms:
        print("\u2714 Chromosomes registered")
    else:
        print("\u2716 Chromosomes not registered")
    if has_seq:
        print("\u2714 Sequence registered")
    else:
        print("\u2716 Sequence not registered")
    if has_tracks:
        print("\u2714 Tracks registered")
    else:
        print("\u2716 Tracks not registered")
    if has_features:
        print("\u2714 Features registered")
    else:
        print("\u2716 Features not registered")

    # print(f"Momics repository  >  {path}")
    print("-".join([""] * 80))
    print(f"{name}")

    if has_chroms:
        print(" \\_ genome")
        print(f"     \\_ chroms     >  {chroms_uri}")

    if has_seq:
        [
            print("     \\_ " + chrom + " ".join([""] * (12 - len(chrom))) + ">  " + mom._build_uri("genome", chrom) + ".tdb")
            for chrom in chrs.iloc[0 : min(2, len(chrs))]
        ]
        if len(chrs) > 2:
            print("     \\_ ...")

    if has_tracks:
        print(" \\_ coverage      ")
        print(f"     \\_ tracks     >  {tracks_uri}")
        [
            print("     \\_ " + chrom + " ".join([""] * (12 - len(chrom))) + ">  " + mom._build_uri("coverage", chrom) + ".tdb")
            for chrom in chrs.iloc[0 : min(2, len(chrs))]
        ]
        if len(chrs) > 2:
            print("     \\_ ...")

    if has_features:
        print(" \\_ features        ")
        print(f"     \\_ features   >  {features_uri}")
        [
            print(
                "     \\_ " + chrom + " ".join([""] * (12 - len(chrom))) + ">  " + mom._build_uri("annotations", chrom) + ".tdb"
            )
            for chrom in chrs.iloc[0 : min(2, len(chrs))]
        ]
        if len(chrs) > 2:
            print("     \\_ ...")

    return True

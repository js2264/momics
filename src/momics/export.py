import logging
import os
from pathlib import Path

import numpy as np
import pyBigWig
import Bio
from Bio import SeqIO

from .momics import Momics
from .multirangequery import MultiRangeQuery
from .utils import _check_feature_name, _check_track_name


def export_track(momics: Momics, track: str, output: Path) -> Momics:
    """Export a track from a `.momics` repository as a `.bw `file.

    Args:
        track (str): Which track to remove
        output (Path): Prefix of the output bigwig file

    Returns:
        Momics: An updated Momics object
    """
    # Abort if `track` is not listed
    _check_track_name(track, momics.tracks())

    # Silence logger
    logging.disable(logging.CRITICAL)

    # Init output file
    bw = pyBigWig.open(output, "w")
    chrom_sizes = momics.chroms()[["chrom", "length"]].apply(tuple, axis=1).tolist()
    bw.addHeader(chrom_sizes)
    for chrom, chrom_length in chrom_sizes:
        q = MultiRangeQuery(momics, chrom).query_tracks(tracks=[track])
        chroms = np.array([chrom] * chrom_length)
        starts = np.array(range(chrom_length))
        ends = starts + 1
        values0 = q.coverage[track][next(iter(q.coverage[track].keys()))]
        bw.addEntries(chroms, starts=starts, ends=ends, values=values0)
    bw.close()


def export_sequence(momics: Momics, output: Path) -> Momics:
    """Export sequence from a `.momics` repository as a `.fa `file.

    Args:
        output (Path): Prefix of the output bigwig file

    Returns:
        Momics: An updated Momics object
    """
    # Silence logger
    logging.disable(logging.CRITICAL)

    if os.path.exists(output):
        os.remove(output)

    # Init output file
    chroms = momics.chroms()["chrom"]
    with open(output, "a") as output_handle:
        for chrom in chroms:
            q = MultiRangeQuery(momics, chrom).query_sequence()
            seq = q.seq["nucleotide"][next(iter(q.seq["nucleotide"].keys()))]
            sr = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id=chrom, description="")
            SeqIO.write(sr, output_handle, "fasta")


def export_features(momics: Momics, features: str, output: Path) -> Momics:
    """Export a features set from a `.momics` repository as a `.bed `file.

    Args:
        features (str): Which features to remove
        output (Path): Prefix of the output BED file

    Returns:
        Momics: An updated Momics object
    """
    # Abort if `features` is not listed
    _check_feature_name(features, momics.features())

    # Init output file
    bed = momics.features(features)
    bed.to_bed(output)
    return True

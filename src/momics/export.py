from pathlib import Path

import numpy as np
import pyBigWig

from .momics import Momics
from .multirangequery import MultiRangeQuery
from .utils import _check_track_name


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

    # Init output file
    bw = pyBigWig.open(output, "w")
    chrom_sizes = momics.chroms()[["chrom", "length"]].apply(tuple, axis=1).tolist()
    bw.addHeader(chrom_sizes)
    for chrom, _ in chrom_sizes:
        q = MultiRangeQuery(momics, chrom).query_tracks().to_df()
        starts = q["position"].to_numpy(dtype=np.int64)
        ends = starts + 1
        values0 = q[track].to_numpy(dtype=np.float32)
        chroms = np.array([chrom] * len(values0))
        bw.addEntries(chroms, starts=starts, ends=ends, values=values0)
    bw.close()

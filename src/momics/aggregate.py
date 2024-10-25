from typing import Literal, Optional
from pathlib import Path
import numpy as np

from .logging import logger
from .utils import dict_to_bigwig


# For this function, the `type` argument can be either "mean" or "sum"
def aggregate(cov, ranges, chrom_sizes, type: Literal["mean", "sum"] = "mean", prefix: Optional[str] = None) -> dict:
    """
    Aggregate query coverage outputs into genome-wide dictionary(ies).
    The coverage over each range is aggregated across all tracks. In the case of
    overlapping ranges, the coverage is averaged.

    Each value of the output dictionary is a dictionary itself, with the keys being the chromosome names
    and the values being the coverage score, averaged for overlapping ranges.

    Args:
        cov (dict): A dictionary of coverage scores, for each track. This is generally the output of
            :func:`MomicsQuery.query_tracks().coverage`.
        ranges (PyRanges): A PyRanges object containing the ranges queried.
        chrom_sizes (dict): A dictionary of chromosome sizes.
        type: The type of aggregation to perform. Can be either "mean" or "sum".
        prefix (str, optional): Prefix to the output `.bw` files to create.
            If provided, queried coverage will be saved for each track in a file
            named `<prefix>_<track_label>.bw`.

    Returns:
        A dictionary of genome-wide coverage scores, for each track. If
        the queried ranges overlap, the coverage is averaged/summed.
        Note that if the output argument is provided, the results for each track will be
        saved to a `<prefix>_<track_label>.bw` file.

    See Also:
        :func:`MomicsQuery.query_tracks()`

    Examples:
        >>> mom = momics.momics.Momics('path/to/momics')
        >>> windows = pr.PyRanges(
        ...     chromosomes = ["I", "I", "I", "I"],
        ...     starts = [0, 5, 10, 20],
        ...     ends = [30, 30, 30, 30],
        ... )
        >>> cov = MomicsQuery(mom, windows).coverage
        >>> aggregate(cov, windows, {"I": 30})
    """
    attrs = cov.keys()
    tracks = {attr: dict() for attr in attrs}

    for attr in iter(attrs):
        attr_cov = cov[attr]
        track = {chrom: np.zeros(size) for chrom, size in chrom_sizes.items()}
        overlap_count = {chrom: np.zeros(size) for chrom, size in chrom_sizes.items()}

        for (_, row), (_, row_cov) in zip(ranges.df.iterrows(), attr_cov.items()):
            chrom = row["Chromosome"]
            start = row["Start"]
            end = row["End"]
            track[chrom][start:end] += row_cov
            overlap_count[chrom][start:end] += 1

        if type == "mean":
            for chrom in track:
                non_zero_mask = overlap_count[chrom] > 0
                track[chrom][non_zero_mask] /= overlap_count[chrom][non_zero_mask]

        tracks[attr] = track

    if prefix is not None:
        bw_paths = []
        for attr in attrs:
            f = Path(f"{prefix}_{attr}.bw")
            p = dict_to_bigwig(tracks[attr], f)
            logger.info(f"Saved coverage for {attr} to {p.name}")
            bw_paths.append(p)
        return tracks

    else:
        return tracks

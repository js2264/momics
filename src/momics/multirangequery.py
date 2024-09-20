import multiprocessing
import os

import numpy as np
import pandas as pd
import tiledb

from .api import Momics
from .utils import parse_ucsc_coordinates


class MultiRangeQuery:
    def __init__(self, momics: Momics, bed: pd.DataFrame):
        """
        Initialize the MultiRangeQuery object.

        Parameters:
        - momics: a Momics object
        - bed: pd.DataFrame with at least three columns `chr`, `start` and `end`.
        """

        self.momics = momics
        if isinstance(bed, str):
            if ":" in bed:
                bed = parse_ucsc_coordinates(bed)
            else:
                chrom = bed
                chroms = self.momics.chroms()
                chrlength = chroms[chroms["chr"] == chrom]["length"][0]
                bed = parse_ucsc_coordinates(f"{chrom}:1-{chrlength}")

        self.bed = bed
        self.coverage = None
        self.seq = None
        groups = self.bed.groupby("chr")
        chrs, indices = np.unique(self.bed["chr"], return_index=True)
        sorted_chrs = chrs[np.argsort(indices)]
        ordered_groups = {key: group for key, group in groups}
        ordered_groups = {key: ordered_groups[key] for key in list(sorted_chrs)}
        all_range_labels = [
            f"{chr}:{start}-{end}"
            for i, (chr, start, end) in enumerate(
                list(zip(self.bed["chr"], self.bed["start"], self.bed["end"]))
            )
        ]
        self.ordered_groups = ordered_groups
        self.all_range_labels = all_range_labels

    def _query_tracks_per_chr(self, chrom, group):

        print(chrom)
        ranges = list(zip(group["start"], group["end"]))

        # Prepare empty long DataFrame without scores, to merge with results
        ranges_str = []
        for _, (start, end) in enumerate(ranges):
            breadth = end - start + 1  # breadth = number of elements in the range
            label = f"{chrom}:{start}-{end}"  # Range1, Range2, Range3, ...
            ranges_str.extend([label] * breadth)  # Repeat the label breadth times
        ranges_df = pd.DataFrame(
            {
                "range": ranges_str,
                "chr": chrom,
                "position": [
                    item for X in ranges for item in np.arange(X[0], X[1] + 1)
                ],
            }
        )

        # Extract scores from tileDB and wrangle them into DataFrame
        ranges_1 = list(zip(group["start"] - 1, group["end"] - 1))
        tdb = os.path.join(self.momics.path, "coverage", f"{chrom}.tdb")
        with tiledb.open(tdb, "r") as A:
            subarray = A.multi_index[ranges_1, :]

        # Reformat to add track and range labels
        tr = self.momics.tracks()
        tr = tr[[x != "None" for x in tr["label"]]]
        subarray_df = pd.merge(tr, pd.DataFrame(subarray), on="idx").drop(
            ["idx", "path"], axis=1
        )
        subarray_df["position"] += 1
        df = pd.merge(ranges_df, subarray_df, on="position", how="left")
        res = {}
        for track in list(tr["label"]):
            res[track] = (
                df[df["label"] == track]
                .groupby("range")["scores"]
                .apply(list)
                .to_dict()
            )
        return res

    def query_tracks(self):
        """
        Query multiple coverage ranges from a Momics repo.

        Returns:
        - A MultiRangeQuery object
        """

        ## Process each chr separately
        args = [(chrom, group) for (chrom, group) in self.ordered_groups.items()]
        multiprocessing.set_start_method("spawn", force=True)
        with multiprocessing.Pool(processes=12) as pool:
            results = pool.starmap(self._query_tracks_per_chr, args)
        res = {}
        tr = self.momics.tracks()
        tr = list(tr[[x != "None" for x in tr["label"]]]["label"])
        for track in tr:
            scores_from_all_chrs = [x[track] for x in results]
            d = {k: v for d in scores_from_all_chrs for k, v in d.items()}
            res[track] = {key: d[key] for key in self.all_range_labels if key in d}

        self.coverage = res
        return self

    def query_sequence(self):
        """
        Query multiple sequence ranges from a Momics repo.

        Returns:
        - A MultiRangeQuery object
        """

        ## Process each chr separately
        args = [(chrom, group) for (chrom, group) in self.ordered_groups.items()]
        multiprocessing.set_start_method("spawn", force=True)
        with multiprocessing.Pool(processes=12) as pool:
            results = pool.starmap(self._query_tracks_per_chr, args)
        res = {}
        tr = self.momics.tracks()
        tr = list(tr[[x != "None" for x in tr["label"]]]["label"])
        for track in tr:
            scores_from_all_chrs = [x[track] for x in results]
            d = {k: v for d in scores_from_all_chrs for k, v in d.items()}
            res[track] = {key: d[key] for key in self.all_range_labels if key in d}

        self.coverage = res
        return self

    def to_df(self):
        # Prepare empty long DataFrame without scores, to merge with results
        cov = self.coverage
        if cov is None:
            raise AttributeError(
                "self.coverage is None. Call `self.query_tracks()` to populate it."
            )

        ranges = list(cov[next(iter(cov.keys()))].keys())
        ranges_str = []
        for _, coords in enumerate(ranges):
            chrom, range_part = coords.split(":")
            start = int(range_part.split("-")[0])
            end = int(range_part.split("-")[1])
            label = [{"chr": chrom, "position": x} for x in range(start, end + 1)]
            ranges_str.extend(label)
        df = pd.DataFrame(ranges_str)

        for track in list(cov.keys()):
            df[track] = [value for sublist in cov[track].values() for value in sublist]
        return df

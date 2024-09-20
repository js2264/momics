import tiledb
import numpy as np
import pandas as pd
import os
import multiprocessing


class MultiChromRangeQuery:
    def __init__(self, path: str, bed: pd.DataFrame, tracks: pd.DataFrame):
        """
        Initialize the MultiChromRangeQuery object.

        Parameters:
        - path: Path to a momics.
        - bed: pd.DataFrame with at least three columns `chr`, `start` and `end`.
        - tracks: pd.DataFrame listing the tracks being recovered.
        """
        self.path = path
        self.bed = bed
        self.tracks = tracks
        groups = bed.groupby("chr")
        chrs, indices = np.unique(bed["chr"], return_index=True)
        sorted_chrs = chrs[np.argsort(indices)]
        ordered_groups = {key: group for key, group in groups}
        ordered_groups = {key: ordered_groups[key] for key in list(sorted_chrs)}
        all_range_labels = [
            f"{chr}:{start}-{end}"
            for i, (chr, start, end) in enumerate(
                list(zip(bed["chr"], bed["start"], bed["end"]))
            )
        ]
        self.ordered_groups = ordered_groups
        self.all_range_labels = all_range_labels

    def _query_tracks_per_chr(self, chrom, group):

        print(chrom)
        ranges = list(zip(group["start"], group["end"]))

        # Prepare empty long DataFrame without scores, to merge with results
        ranges_str = []
        for i, (start, end) in enumerate(ranges):
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
        tdb = os.path.join(self.path, "coverage", f"{chrom}.tdb")
        with tiledb.open(tdb, "r") as A:
            subarray = A.multi_index[ranges_1, :]

        # Reformat to add track and range labels
        subarray_df = pd.merge(self.tracks, pd.DataFrame(subarray), on="idx").drop(
            ["idx"], axis=1
        )
        subarray_df["position"] += 1
        df = pd.merge(ranges_df, subarray_df, on="position", how="left")
        l = {}
        for track in list(self.tracks["label"]):
            l[track] = (
                df[df["label"] == track]
                .groupby("range")["scores"]
                .apply(list)
                .to_dict()
            )
        return l

    def query_tracks(self):
        """
        Query multiple coverage ranges from a Momics repo.

        Returns:
        - A DataFrame
        """

        ## Process each chr separately
        args = [(chrom, group) for (chrom, group) in self.ordered_groups.items()]
        multiprocessing.set_start_method("spawn", force=True)
        with multiprocessing.Pool(processes=12) as pool:
            results = pool.starmap(self._query_tracks_per_chr, args)

        l = {}
        for track in list(self.tracks["label"]):
            scores_from_all_chrs = [x[track] for x in results]
            d = {k: v for d in scores_from_all_chrs for k, v in d.items()}
            l[track] = {key: d[key] for key in self.all_range_labels if key in d}

        return l

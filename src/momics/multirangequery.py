import concurrent.futures
from pathlib import Path

import numpy as np
import pandas as pd
import json
import pickle
import tiledb
import threading
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .momics import Momics
from .utils import parse_ucsc_coordinates
from .logging import logger

lock = threading.Lock()


class MultiRangeQuery:
    """A class to query `.momics` repositories.

    Attributes
    ----------
    momics (Momics): a local `.momics` repository.
    queries (dict): Dict. of pd.DataFrames with at least three columns `chrom`, `start` and `end`, one per chromosome.
    query_labels (list): List of UCSC-style coordinates.
    coverage (dict): Dictionary of coverage scores extracted from the `.momics` repository, populated after calling `q.query_tracks()`
    seq (dict): Dictionary of sequences extracted from the `.momics` repository, populated after calling `q.query_seq()`
    """

    def __init__(self, momics: Momics, bed: pd.DataFrame):
        """Initialize the MultiRangeQuery object.

        Args:
            momics (Momics): a Momics object
            bed (pd.DataFrame): pd.DataFrame with at least three columns `chrom`, `start` and `end`.
        """
        self.momics = momics
        if isinstance(bed, str):
            if ":" in bed:
                bed = parse_ucsc_coordinates(bed)
            else:
                chrom = bed
                chroms = self.momics.chroms()
                chrlength = chroms[chroms["chrom"] == chrom]["length"].iloc[0]
                bed = parse_ucsc_coordinates(f"{chrom}:1-{chrlength}")

        groups = bed.groupby("chrom")
        chrs, indices = np.unique(bed["chrom"], return_index=True)
        sorted_chrs = chrs[np.argsort(indices)]
        queries = {key: group for key, group in groups}
        queries = {key: queries[key] for key in list(sorted_chrs)}
        self.queries = queries
        query_labels = [
            f"{chr}:{start}-{end}"
            for _, (chr, start, end) in enumerate(
                list(zip(bed["chrom"], bed["start"], bed["end"]))
            )
        ]
        self.query_labels = query_labels
        self.coverage = None
        self.seq = None

    def query_tracks(self, threads: int = 1) -> "MultiRangeQuery":
        """Query multiple coverage ranges from a Momics repo.

        Args:
            threads (int, optional): Number of threads for parallel query. Defaults to 1.

        Returns:
            MultiRangeQuery: MultiRangeQuery: An updated MultiRangeQuery object
        """

        def _query_tracks_per_chr(self, chrom, group, tracks):

            ranges = list(zip(group["start"], group["end"]))

            # Prepare empty long DataFrame without scores, to merge with results
            ranges_str = []
            for _, (start, end) in enumerate(ranges):
                breadth = end - start + 1
                label = f"{chrom}:{start}-{end}"
                ranges_str.extend([label] * breadth)
            ranges_df = pd.DataFrame(
                {
                    "range": ranges_str,
                    "chrom": chrom,
                    "position": [
                        item for X in ranges for item in np.arange(X[0], X[1] + 1)
                    ],
                }
            )

            # Extract scores from tileDB and wrangle them into DataFrame
            ranges_1 = list(zip(group["start"] - 1, group["end"] - 1))
            tdb = self.momics._build_uri("coverage", f"{chrom}.tdb")
            with tiledb.open(tdb, "r", ctx=self.momics.cfg.ctx) as A:
                subarray = A.multi_index[ranges_1, :]

            # Reformat to add track and range labels
            tr = tracks[[x != "None" for x in tracks["label"]]]
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

        def _log_task_completion(future, chrom, ntasks, completed_tasks):
            if future.exception() is not None:
                logger.error(
                    f"Querying tracks for chromosome {chrom} failed with exception: {future.exception()}"
                )
            else:
                with lock:
                    completed_tasks[0] += 1
                logger.info(
                    f"task {completed_tasks[0]}/{ntasks} :: Queried tracks for chromosome {chrom}."
                )

        tracks = self.momics.tracks()
        tasks = [(chrom, group) for (chrom, group) in self.queries.items()]
        ntasks = len(tasks)
        completed_tasks = [0]
        threads = min(threads, ntasks)

        if threads == 1:
            results = []
            for chrom, group in tasks:
                results.append(_query_tracks_per_chr(self, chrom, group, tracks))
                completed_tasks[0] += 1
                logger.info(
                    f"task {completed_tasks[0]}/{ntasks} :: Queried tracks for chromosome {chrom}."
                )

        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
                futures = []
                for chrom, group in tasks:
                    future = executor.submit(
                        _query_tracks_per_chr, self, chrom, group, tracks
                    )
                    future.add_done_callback(
                        lambda f, c=chrom: _log_task_completion(
                            f, c, ntasks, completed_tasks
                        )
                    )
                    futures.append(future)
                concurrent.futures.wait(futures)

            results = []
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())

        res = {}
        tr = list(tracks[[x != "None" for x in tracks["label"]]]["label"])
        for track in tr:
            scores_from_all_chrs = [x[track] for x in results]
            d = {k: v for d in scores_from_all_chrs for k, v in d.items()}
            res[track] = {key: d[key] for key in self.query_labels if key in d}

        self.coverage = res
        return self

    def query_sequence(self, threads: int = 1) -> "MultiRangeQuery":
        """Query multiple sequence ranges from a Momics repo.

        Args:
            threads (int, optional): Number of threads for parallel query. Defaults to 1.

        Returns:
            MultiRangeQuery: An updated MultiRangeQuery object
        """

        def _query_seq_per_chr(self, chrom, group):
            ranges = list(zip(group["start"], group["end"]))
            seqs = {}
            tdb = self.momics._build_uri("genome", "sequence", f"{chrom}.tdb")
            for _, (start, end) in enumerate(ranges):
                with tiledb.open(tdb, "r", ctx=self.momics.cfg.ctx) as A:
                    seq = A.df[(start - 1) : (end - 1)]["nucleotide"]
                seqs[f"{chrom}:{start}-{end}"] = "".join(seq)
            return seqs

        def _log_task_completion(future, chrom, ntasks, completed_tasks):
            if future.exception() is not None:
                logger.error(
                    f"Querying sequences for chromosome {chrom} failed with exception: {future.exception()}"
                )
            else:
                with lock:
                    completed_tasks[0] += 1
                logger.info(
                    f"task {completed_tasks[0]}/{ntasks} :: Queried sequences for chromosome {chrom}."
                )

        tasks = [(chrom, group) for (chrom, group) in self.queries.items()]
        ntasks = len(tasks)
        completed_tasks = [0]
        threads = min(threads, ntasks)
        if threads == 1:
            seqs = []
            for chrom, group in tasks:
                seqs.append(_query_seq_per_chr(self, chrom, group))
                completed_tasks[0] += 1
                logger.info(
                    f"task {completed_tasks[0]}/{ntasks} :: Queried sequences for chromosome {chrom}."
                )
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
                futures = []
                for chrom, group in tasks:
                    future = executor.submit(_query_seq_per_chr, self, chrom, group)
                    future.add_done_callback(
                        lambda f, c=chrom: _log_task_completion(
                            f, c, ntasks, completed_tasks
                        )
                    )
                    futures.append(future)
                seqs = []
                for future in concurrent.futures.as_completed(futures):
                    seqs.append(future.result())

        mseqs = {}
        for d in seqs:
            mseqs.update(d)
        self.seq = {"seq": mseqs}
        return self

    def to_df(self) -> pd.DataFrame:
        """Parse self.coverage attribute to a pd.DataFrame

        Returns:
            pd.DataFrame: `self.coverage` dictionary wrangled into a pd.DataFrame
        """
        # Prepare empty long DataFrame without scores, to merge with results
        cov = self.coverage
        if cov is None:
            raise AttributeError(
                "self.coverage is None. Call `self.query_tracks()` to populate it."
            )

        ranges_str = []
        for _, coords in enumerate(self.query_labels):
            chrom, range_part = coords.split(":")
            start = int(range_part.split("-")[0])
            end = int(range_part.split("-")[1])
            label = [
                {"range": coords, "chrom": chrom, "position": x}
                for x in range(start, end + 1)
            ]
            ranges_str.extend(label)
        df = pd.DataFrame(ranges_str)

        for track in list(cov.keys()):
            df[track] = [value for sublist in cov[track].values() for value in sublist]
        return df

    def to_fa(self) -> SeqRecord:
        """Parse self.seq attribute to a SeqRecord

        Returns:
            SeqRecord: `self.seq` dictionary wrangled into a SeqRecord
        """

        seq = self.seq
        if seq is None:
            raise AttributeError(
                "self.seq is None. Call `self.query_sequence()` to populate it."
            )

        seq_records = []
        for header, sequence in seq["seq"].items():
            seq_record = SeqRecord(Seq(sequence), id=header, description="")
            seq_records.append(seq_record)

        return seq_records

    def to_npz(self, output: Path):
        """Write the results of a multi-range query to a NPZ file.

        Args:
            output (Path): Path to the output NPZ file.
        """
        if self.coverage is None:
            raise AttributeError(
                "self.coverage is None. Call `self.query_tracks()` to populate it."
            )
        if self.seq is None:
            raise AttributeError(
                "self.seq is None. Call `self.query_sequence()` to populate it."
            )
        serialized_cov = pickle.dumps(self.coverage)
        serialized_seq = pickle.dumps(self.seq)
        logger.info(f"Saving results of multi-range query to {output}...")
        with open(output, "wb") as f:
            np.savez_compressed(f, coverage=serialized_cov, seq=serialized_seq)

    def to_json(self, output: Path):
        """Write the results of a multi-range query to a JSON file.

        Args:
            output (Path): Path to the output JSON file.
        """
        data = self.coverage
        data["seq"] = self.seq["seq"]
        logger.info(f"Saving results of multi-range query to {output}...")
        with open(output, "w") as json_file:
            json.dump(data, json_file, indent=4)

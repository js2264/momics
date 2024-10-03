import collections
import concurrent.futures
import os
import time
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
    coordinates (list): List of UCSC-style coordinates.
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

        ranges = [
            f"{chr}:{start}-{end}"
            for _, (chr, start, end) in enumerate(
                list(zip(bed["chrom"], bed["start"], bed["end"]))
            )
        ]
        self.ranges = ranges
        self.coverage = None
        self.seq = None

    @staticmethod
    def _query_tracks_per_chr(base_uri, ranges, cfg_dict):

        try:

            start0 = time.time()
            ranges_per_chr = collections.defaultdict(list)
            for r in ranges:
                chrom, value = r.split(":")
                ranges_per_chr[chrom].append(value)
            logger.info(f"Generating dict per chroms :: {time.time() - start0}")

            start0 = time.time()
            scores = collections.defaultdict(list)
            for chrom, subranges in ranges_per_chr.items():
                scores[chrom] = collections.defaultdict(list)
                tdb = os.path.join(base_uri, "coverage", f"{chrom}.tdb")

                # Extract scores from tileDB and wrangle them into DataFrame
                query = [
                    slice(int(start), int(end) - 1)
                    for start, end in (r.split("-") for r in subranges)
                ]
                with tiledb.open(tdb, "r", config=tiledb.Config(cfg_dict)) as A:
                    subarray = A.multi_index[query,]

                # Extract scores from tileDB and wrangle them into DataFrame
                attrs = list(subarray.keys())
                for attr in attrs:
                    cov = subarray[attr]
                    coverage_per_slice = []
                    start_idx = 0
                    query_lengths = [s.stop - s.start for s in query]
                    for length in query_lengths:
                        coverage_per_slice.append(cov[start_idx : start_idx + length])
                        start_idx += length
                    scores[chrom][attr].append(coverage_per_slice)

            logger.info(f"Extracting scores :: {time.time() - start0}")

            # `scores` are dictionaries of keys: chroms, values: list of attributes
            # Each subdirectory is a dict of keys: attributes, values: list of coverage scores
            starto = time.time()
            combined_scores = collections.defaultdict(list)
            for chrom, attrs in scores.items():
                for attr, scores_list in attrs.items():
                    combined_scores[attr].extend(scores_list[0])
            logger.info(f"Reformatting scores :: {time.time() - start0}")

            return combined_scores

        except Exception as e:
            logger.error(f"Error processing query chunk: {e}")
            raise

    def query_tracks(self, threads: int = 1) -> "MultiRangeQuery":
        """Query multiple coverage ranges from a Momics repo.

        Args:
            threads (int, optional): Number of threads for parallel query. Defaults to 1.

        Returns:
            MultiRangeQuery: MultiRangeQuery: An updated MultiRangeQuery object
        """

        def _log_task_completion(future, ntasks, completed_tasks):
            if future.exception() is not None:
                logger.error(
                    f"Querying tracks [task {completed_tasks+1}] failed with exception: {future.exception()}"
                )
            else:
                with lock:
                    completed_tasks[0] += 1
                logger.info(f"task {completed_tasks[0]}/{ntasks} finished.")

        cfg_dict = {k: v for k, v in self.momics.cfg.cfg.items()}
        tasks = [
            self.ranges[
                i * len(self.ranges) // threads : (i + 1) * len(self.ranges) // threads
            ]
            for i in range(threads)
        ]
        tasks = [t for t in tasks if len(t) > 0]
        ntasks = len(tasks)
        completed_tasks = [0]
        threads = min(threads, ntasks)
        start0 = time.time()

        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            futures = []
            for r in tasks:
                future = executor.submit(
                    self._query_tracks_per_chr, self.momics.path, r, cfg_dict
                )
                future.add_done_callback(
                    lambda f: _log_task_completion(f, ntasks, completed_tasks)
                )
                futures.append(future)
            concurrent.futures.wait(futures)

        results = []
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())

        combined_results = collections.defaultdict(dict)
        attrs = list(results[0].keys())
        for attr in attrs:
            combined_results[attr] = []
        for i in range(len(results)):
            for attr in attrs:
                combined_results[attr].extend(results[i][attr])

        self.coverage = combined_results
        t = time.time() - start0
        logger.info(f"Query completed in {round(t,4)}s.")
        return self

    @staticmethod
    def _query_seq_per_chr(chrom, group, tdb, cfg_dict):
        start0 = time.time()
        ranges = list(zip(group["start"], group["end"]))
        seqs = {}
        with tiledb.open(tdb, "r", ctx=tiledb.Ctx(tiledb.Config(cfg_dict))) as A:
            for _, (start, end) in enumerate(ranges):
                seq = A.df[(start - 1) : (end - 1)]["nucleotide"]
                seqs[f"{chrom}:{start}-{end}"] = "".join(seq)
        logger.debug(f"Fetching sequence :: {time.time() - start0}")
        return seqs

    def query_sequence(self, threads: int = 1) -> "MultiRangeQuery":
        """Query multiple sequence ranges from a Momics repo.

        Args:
            threads (int, optional): Number of threads for parallel query. Defaults to 1.

        Returns:
            MultiRangeQuery: An updated MultiRangeQuery object
        """

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

        tasks = [
            (
                chrom,
                group,
                self.momics._build_uri("genome", "sequence", f"{chrom}.tdb"),
                {k: v for k, v in self.momics.cfg.cfg.items()},
            )
            for (chrom, group) in self.queries.items()
        ]
        ntasks = len(tasks)
        completed_tasks = [0]
        threads = min(threads, ntasks)
        start0 = time.time()

        if threads == 1:
            seqs = []
            for chrom, group, tdb, cfg_dict in tasks:
                seqs.append(self._query_seq_per_chr(chrom, group, tdb, cfg_dict))
                completed_tasks[0] += 1
                logger.info(
                    f"task {completed_tasks[0]}/{ntasks} :: Queried sequences for chromosome {chrom}."
                )
        else:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=threads
            ) as executor:
                futures = []
                for chrom, group, tdb, cfg_dict in tasks:
                    future = executor.submit(
                        self._query_seq_per_chr, chrom, group, tdb, cfg_dict
                    )
                    future.add_done_callback(
                        lambda f, c=chrom: _log_task_completion(
                            f, c, ntasks, completed_tasks
                        )
                    )
                    futures.append(future)
                concurrent.futures.wait(futures)

            seqs = []
            for future in concurrent.futures.as_completed(futures):
                seqs.append(future.result())

        mseqs = {}
        for d in seqs:
            mseqs.update(d)
        self.seq = {"seq": mseqs}
        t = time.time() - start0
        logger.info(f"Query completed in {round(t,4)}s.")
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
        for _, coords in enumerate(self.coordinates):
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

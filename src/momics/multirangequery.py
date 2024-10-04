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
        if not isinstance(momics, Momics):
            raise ValueError("momics must be a `Momics` object.")
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
    def _query_tracks_per_batch(base_uri, ranges, cfg_dict, tracks):
        try:
            # Split ranges by chromosome
            ranges_per_chr = collections.defaultdict(list)
            for r in ranges:
                chrom, value = r.split(":")
                ranges_per_chr[chrom].append(value)

            # Get attributes
            _c = list(ranges_per_chr.keys())[0]
            _sch = tiledb.open(
                os.path.join(base_uri, "coverage", f"{_c}.tdb"),
                "r",
                config=tiledb.Config(cfg_dict),
            ).schema
            attrs = [_sch.attr(i).name for i in range(_sch.nattr)]
            if tracks is not None:
                attrs = [attr for attr in attrs if attr in tracks]

            # Prepare empty dictionary {attr1: { ranges1: ..., ranges2: ... }, attr2: {}, ...}
            results = {attr: collections.defaultdict(list) for attr in attrs}

            for chrom, subranges in ranges_per_chr.items():
                keys = [f"{chrom}:{i}" for i in subranges]
                tdb = os.path.join(base_uri, "coverage", f"{chrom}.tdb")

                # Extract scores from tileDB and wrangle them into DataFrame
                query = [
                    slice(int(start) - 1, int(end))
                    for start, end in (r.split("-") for r in subranges)
                ]
                cfg = tiledb.Config(cfg_dict)
                cfg.update({"sm.compute_concurrency_level": 1})
                cfg.update({"sm.io_concurrency_level": 1})
                with tiledb.open(tdb, "r", config=cfg) as A:
                    # subarray = A.multi_index[query,]
                    subarray = A.query(attrs=attrs).multi_index[query,]

                # Extract scores from tileDB and wrangle them into DataFrame
                # This is the tricky bit, because tileDB returns a dict of attributes
                # and for each attribute, there is only a single list of scores
                # all concatenated together. We need to split them back into the
                # original slices.
                for attr in attrs:
                    cov = subarray[attr]
                    start_idx = 0
                    query_lengths = [s.stop - s.start for s in query]
                    for i, length in enumerate(query_lengths):
                        results[attr][keys[i]] = cov[start_idx : start_idx + length]
                        start_idx += length

            return results

        except Exception as e:
            logger.error(f"Error processing query batch: {e}")
            raise

    def query_tracks(self, threads: int = 1, tracks: list = None) -> "MultiRangeQuery":
        """Query multiple coverage ranges from a Momics repo.

        Args:
            threads (int, optional): Number of threads for parallel query. Defaults to 1.
            tracks (list, optional): List of tracks to query. Defaults to None, which queries all tracks.

        Returns:
            MultiRangeQuery: MultiRangeQuery: An updated MultiRangeQuery object
        """

        def _log_task_completion(future, ntasks, completed_tasks):
            if future.exception() is not None:
                logger.error(
                    f"Querying tracks [task {completed_tasks[0]+1}] failed with exception: {future.exception()}"
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

        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            futures = []
            for i, r in enumerate(tasks):
                future = executor.submit(
                    self._query_tracks_per_batch, self.momics.path, r, cfg_dict, tracks
                )
                future.add_done_callback(
                    lambda f: _log_task_completion(f, ntasks, completed_tasks)
                )
                futures.append((i, future))
            concurrent.futures.wait([f[1] for f in futures])

        results = [None] * len(tasks)
        for i, future in futures:
            results[i] = future.result()

        attrs = results[0].keys()
        combined_results = {attr: {} for attr in attrs}
        for d in results:
            for attr in attrs:
                combined_results[attr].update(d[attr])

        self.coverage = combined_results
        t = time.time() - start0
        logger.info(f"Query completed in {round(t,4)}s.")
        return self

    @staticmethod
    def _query_sequence_per_batch(base_uri, ranges, cfg_dict):
        try:
            start0 = time.time()
            ranges_per_chr = collections.defaultdict(list)
            for r in ranges:
                chrom, value = r.split(":")
                ranges_per_chr[chrom].append(value)
            logger.debug(f"Generating dict per chroms :: {time.time() - start0}")

            start0 = time.time()
            seqs = collections.defaultdict(list)
            for chrom, subranges in ranges_per_chr.items():
                keys = [f"{chrom}:{i}" for i in subranges]
                seqs[chrom] = collections.defaultdict(list)
                tdb = os.path.join(base_uri, "genome", "sequence", f"{chrom}.tdb")

                # Extract sequence from tileDB
                query = [
                    slice(int(start) - 1, int(end))
                    for start, end in (r.split("-") for r in subranges)
                ]
                with tiledb.open(tdb, "r", config=tiledb.Config(cfg_dict)) as A:
                    subarray = A.multi_index[query,]

                # Extract seqs from tileDB and wrangle them into DataFrame
                # This is the tricky bit, because tileDB returns a dict of attributes
                # and for each attribute, there is only a single list of seqs
                # all concatenated together. We need to split them back into the
                # original slices.
                seq = subarray["nucleotide"]
                seq_per_slice = []
                start_idx = 0
                query_lengths = [s.stop - s.start for s in query]
                for length in query_lengths:
                    x = seq[start_idx : start_idx + length]
                    seq_per_slice.append("".join([nt.decode("utf-8") for nt in x]))
                    start_idx += length
                seqs[chrom] = {keys[i]: seq_per_slice[i] for i in range(len(keys))}

            logger.debug(f"Extracting sequences :: {time.time() - start0}")

            return dict(seqs)

        except Exception as e:
            logger.error(f"Error processing query batch: {e}")
            raise

    def query_sequence(self, threads: int = 1) -> "MultiRangeQuery":
        """Query multiple sequence ranges from a Momics repo.

        Args:
            threads (int, optional): Number of threads for parallel query. Defaults to 1.

        Returns:
            MultiRangeQuery: An updated MultiRangeQuery object
        """

        def _log_task_completion(future, ntasks, completed_tasks):
            if future.exception() is not None:
                logger.error(
                    f"Querying sequence [task {completed_tasks[0]+1}] failed with exception: {future.exception()}"
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
            for i, r in enumerate(tasks):
                future = executor.submit(
                    self._query_sequence_per_batch, self.momics.path, r, cfg_dict
                )
                future.add_done_callback(
                    lambda f: _log_task_completion(f, ntasks, completed_tasks)
                )
                futures.append((i, future))
            concurrent.futures.wait([f[1] for f in futures])

        results = [None] * len(tasks)
        for i, future in futures:
            results[i] = future.result()

        combined_results = {}
        for d in results:
            for chrom_dict in d.values():
                combined_results.update(chrom_dict)

        self.seq = {"nucleotide": combined_results}
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
        for _, coords in enumerate(self.ranges):
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
        for header, sequence in seq["nucleotide"].items():
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
        for key, _ in data.items():
            for key2, value2 in data[key].items():
                data[key][key2] = value2.tolist()
        data["nucleotide"] = self.seq["nucleotide"]
        logger.info(f"Saving results of multi-range query to {output}...")
        with open(output, "w") as json_file:
            json.dump(data, json_file, indent=4)

    def dump(self, output: Path):
        """Write the results of a multi-range query to a JSON file.

        Args:
            output (Path): Path to the output JSON file.
        """
        data = self.coverage
        for key, _ in data.items():
            for key2, value2 in data[key].items():
                data[key][key2] = value2.tolist()
        data["nucleotide"] = self.seq["nucleotide"]
        logger.info(f"Saving results of multi-range query to {output}...")
        with open(output, "w") as json_file:
            json.dump(data, json_file, indent=4)

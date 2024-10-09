import collections
import json
import pickle
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import psutil
import pybedtools
import tiledb
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .logging import logger
from .momics import Momics
from .utils import parse_ucsc_coordinates


class MultiRangeQuery:
    """A class to query `.momics` repositories.

    Attributes
    ----------
    momics (Momics): a local `.momics` repository.
    queries (dict): Dict. of pd.DataFrames with at least three columns \
        `chrom`, `start` and `end`, one per chromosome.
    coordinates (list): List of UCSC-style coordinates.
    coverage (dict): Dictionary of coverage scores extracted from the \
        `.momics` repository, populated after calling `q.query_tracks()`
    seq (dict): Dictionary of sequences extracted from the `.momics` \
        repository, populated after calling `q.query_seq()`
    """

    def __init__(self, momics: Momics, bed: pybedtools.BedTool):
        """Initialize the MultiRangeQuery object.

        Args:
            momics (Momics): a Momics object
            bed (pd.DataFrame): pybedtools.BedTool object
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

        # ranges = [f"{inter.chrom}:{inter.start}-{inter.end}" for inter in bed]
        self.ranges = bed
        self.coverage = None
        self.seq = None

    def _check_memory_available(self, n):
        estimated_required_memory = (4 * n * sum([s.stop - s.start + 1 for s in self.ranges])) * 1.2
        emem = round(estimated_required_memory / 1e9, 2)
        avail_mem = psutil.virtual_memory().available
        amem = round(avail_mem / 1e9, 2)
        if estimated_required_memory > avail_mem:
            logger.warning(
                f"Estimated required memory ({emem}GB) exceeds available memory \
                    ({amem}GB)."
            )

    def _query_tracks_per_batch(self, chrom, ranges, attrs, cfg):
        try:

            # Prepare queries: list of slices [(start, stop), (start, stop), ...]
            start0 = time.time()
            query = [slice(int(start) - 1, int(end)) for start, end in (r.split("-") for r in ranges)]
            logger.debug(f"define query in {round(time.time() - start0,4)}s")

            # Query tiledb
            start0 = time.time()
            tdb = self.momics._build_uri("coverage", f"{chrom}.tdb")
            with tiledb.open(tdb, "r", config=cfg) as A:
                subarray = A.query(attrs=attrs).multi_index[
                    query,
                ]
            logger.debug(f"query tiledb in {round(time.time() - start0,4)}s")

            # Extract scores from tileDB and wrangle them into DataFrame
            # This is the tricky bit, because tileDB returns a dict of attributes
            # and for each attribute, there is only a single list of scores
            # all concatenated together. We need to split them back into the
            # original slices.
            start0 = time.time()
            results = {attr: collections.defaultdict(list) for attr in attrs}
            keys = [f"{chrom}:{i}" for i in ranges]
            for attr in attrs:
                cov = subarray[attr]
                start_idx = 0
                query_lengths = [s.stop - s.start for s in query]
                for i, length in enumerate(query_lengths):
                    results[attr][keys[i]] = cov[start_idx : start_idx + length]
                    start_idx += length
            logger.debug(f"wrangle data in {round(time.time() - start0,4)}s")

            return results

        except Exception as e:
            logger.error(f"Error processing query batch: {e}")
            raise

    def query_tracks(self, threads: Optional[int] = None, tracks: Optional[list] = None) -> "MultiRangeQuery":
        """Query multiple coverage ranges from a Momics repo.

        Args:
            threads (int, optional): Number of threads for parallel query. \
                Defaults to all.
            tracks (list, optional): List of tracks to query. Defaults to None, \
                which queries all tracks.

        Returns:
            MultiRangeQuery: MultiRangeQuery: An updated MultiRangeQuery object
        """

        start0 = time.time()

        # Limit tiledb threads
        cfg = self.momics.cfg.cfg
        if threads is not None:
            cfg.update({"sm.compute_concurrency_level": threads})
            cfg.update({"sm.io_concurrency_level": threads})

        # Extract attributes from schema
        chroms = list(dict.fromkeys([x.chrom for x in self.ranges]))
        _sch = tiledb.open(
            self.momics._build_uri("coverage", f"{chroms[0]}.tdb"),
            "r",
            config=cfg,
        ).schema
        attrs = [_sch.attr(i).name for i in range(_sch.nattr)]
        if tracks is not None:
            attrs = [attr for attr in attrs if attr in tracks]

        # Check memory available and warn if it's not enough
        self._check_memory_available(len(attrs))

        # Split ranges by chromosome
        ranges_per_chrom = collections.defaultdict(list)
        for r in self.ranges:
            chrom = r.chrom
            ranges_per_chrom[chrom].append(f"{r.start}-{r.end}")

        # Prepare empty dictionary of results {attr1: { ranges1: ., ranges2: .}, ...}
        results = []
        for chrom in chroms:
            logger.debug(chrom)
            results.append(
                self._query_tracks_per_batch(
                    chrom=chrom,
                    ranges=ranges_per_chrom[chrom],
                    attrs=attrs,
                    cfg=cfg,
                )
            )

        combined_results = {attr: {} for attr in attrs}
        for d in results:
            for attr in attrs:
                combined_results[attr].update(d[attr])

        self.coverage = combined_results
        t = time.time() - start0
        logger.info(f"Query completed in {round(t,4)}s.")
        return self

    def _query_seq_per_batch(self, chrom, ranges, attrs, cfg):
        try:

            # Prepare queries: list of slices [(start, stop), (start, stop), ...]
            start0 = time.time()
            query = [slice(int(start) - 1, int(end)) for start, end in (r.split("-") for r in ranges)]
            logger.debug(f"define query in {round(time.time() - start0,4)}s")

            # Query tiledb
            start0 = time.time()
            tdb = self.momics._build_uri("genome", "sequence", f"{chrom}.tdb")
            with tiledb.open(tdb, "r", config=cfg) as A:
                subarray = A.multi_index[
                    query,
                ]
            logger.debug(f"query tiledb in {round(time.time() - start0,4)}s")

            # Extract scores from tileDB and wrangle them into DataFrame
            # This is the tricky bit, because tileDB returns a dict of attributes
            # and for each attribute, there is only a single list of scores
            # all concatenated together. We need to split them back into the
            # original slices.
            start0 = time.time()
            results = {attr: collections.defaultdict(list) for attr in attrs}
            keys = [f"{chrom}:{i}" for i in ranges]
            for attr in attrs:
                seq = subarray[attr]
                start_idx = 0
                query_lengths = [s.stop - s.start for s in query]
                for i, length in enumerate(query_lengths):
                    results[attr][keys[i]] = "".join(seq[start_idx : start_idx + length].tolist())
                    start_idx += length
            logger.debug(f"wrangle data in {round(time.time() - start0,4)}s")

            return dict(results)

        except Exception as e:
            logger.error(f"Error processing query batch: {e}")
            raise

    def query_sequence(self, threads: Optional[int] = None) -> "MultiRangeQuery":
        """Query multiple sequence ranges from a Momics repo.

        Args:
            threads (int, optional): Number of threads for parallel query. \
                Defaults to all.

        Returns:
            MultiRangeQuery: An updated MultiRangeQuery object
        """

        start0 = time.time()

        # Limit tiledb threads
        cfg = self.momics.cfg.cfg
        if threads is not None:
            cfg.update({"sm.compute_concurrency_level": threads})
            cfg.update({"sm.io_concurrency_level": threads})

        # Split ranges by chromosome
        attrs = ["nucleotide"]
        ranges_per_chrom = collections.defaultdict(list)
        for r in self.ranges:
            chrom = r.chrom
            ranges_per_chrom[chrom].append(f"{r.start}-{r.end}")

        # Prepare empty dictionary of results {attr1: { ranges1: ., ranges2: . }, .}
        results = []
        for chrom in ranges_per_chrom.keys():
            logger.debug(chrom)
            results.append(
                self._query_seq_per_batch(
                    chrom=chrom,
                    ranges=ranges_per_chrom[chrom],
                    attrs=attrs,
                    cfg=cfg,
                )
            )

        combined_results = {attr: {} for attr in attrs}
        for d in results:
            for attr in attrs:
                combined_results[attr].update(d[attr])

        self.seq = combined_results
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
            raise AttributeError("self.coverage is None. Call `self.query_tracks()` to populate it.")

        ranges_str = []
        for inter in self.ranges:
            chrom = inter.chrom
            start = inter.start
            end = inter.end
            label = [{"range": f"{chrom}:{start}-{end}", "chrom": chrom, "position": x} for x in range(start, end + 1)]
            ranges_str.extend(label)
        df = pd.DataFrame(ranges_str)

        for track in list(cov.keys()):
            df[track] = [value for sublist in cov[track].values() for value in sublist]
        return df

    def to_SeqRecord(self) -> SeqRecord:
        """Parse self.seq attribute to a SeqRecord

        Returns:
            SeqRecord: `self.seq` dictionary wrangled into a SeqRecord
        """

        seq = self.seq
        if seq is None:
            raise AttributeError("self.seq is None. Call `self.query_sequence()` to populate it.")

        seq_records = []
        for header, sequence in seq["nucleotide"].items():
            # seq_string = "".join(sequence.astype(str))
            seq_string = "".join(sequence)
            seq_record = SeqRecord(Seq(seq_string), id=header, description="")
            seq_records.append(seq_record)

        return seq_records

    def to_npz(self, output: Path):
        """Write the results of a multi-range query to a NPZ file.

        Args:
            output (Path): Path to the output NPZ file.
        """
        if self.coverage is None:
            raise AttributeError("self.coverage is None. Call `self.query_tracks()` to populate it.")
        if self.seq is None:
            raise AttributeError("self.seq is None. Call `self.query_sequence()` to populate it.")
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

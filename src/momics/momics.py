from datetime import datetime
import os
import tempfile
from time import sleep
import time
from typing import Dict, Optional
from pathlib import Path
import concurrent.futures
import threading

import numpy as np
import pandas as pd
import pyBigWig
import pyfaidx
import tiledb

from . import utils
from .config import MomicsConfig
from .logging import logger


lock = threading.Lock()


class Momics:
    """
    A class to manipulate `.momics` repositories.

    Attributes
    ----------
    path : str
        Path to a `.momics` repository.
    """

    def __init__(
        self,
        path: str,
        config: Optional[MomicsConfig] = None,
    ):
        """
        Initialize the Momics class.
        By default, a `.momics` repository is created at the specified path if it does not already exist.

        Args:
            path (str): Path to a `.momics` repository.
        """

        self.path = path
        if config is None:
            config = MomicsConfig()
        self.cfg = config

        ## Check if folder exists. If not, create it.
        if not self.cfg.vfs.is_dir(self.path):
            self.cfg.vfs.create_dir(self.path)
            self._create_repository()
            logger.info(f"Created {self.path}")
        else:
            logger.info(f"Found {self.path}")

    def _is_cloud_hosted(self):
        if self.path.startswith(("s3://", "gcs://", "azure://")):
            return self.path.split("://")[0]
        else:
            return False

    def _build_uri(self, *subdirs: str) -> str:
        if self._is_cloud_hosted():
            return "/".join([self.path.rstrip("/")] + list(subdirs))
        else:
            return str(Path(self.path).joinpath(*subdirs))

    def _create_repository(self):
        genome_path = self._build_uri("genome")
        seq_path = self._build_uri("genome", "sequence")
        coverage_path = self._build_uri("coverage")
        features_path = self._build_uri("features")
        tiledb.group_create(self.path, ctx=self.cfg.ctx)
        tiledb.group_create(genome_path, ctx=self.cfg.ctx)
        tiledb.group_create(coverage_path, ctx=self.cfg.ctx)
        tiledb.group_create(features_path, ctx=self.cfg.ctx)
        tiledb.group_create(seq_path, ctx=self.cfg.ctx)

    def _get_table(self, uri: str) -> Optional[pd.DataFrame]:
        if not self.cfg.vfs.is_dir(uri):
            raise FileExistsError(f"{uri} does not exist.")

        with tiledb.open(uri, "r", ctx=self.cfg.ctx) as A:
            if A.schema.sparse:
                a = A.df[:]
            else:
                a = A.df[0 : len(A) - 1]

        return a

    def _create_sequence_schema(self, tile: int, compression: int):
        # Create every /sequence/{chrom}.tdb
        chroms = self.chroms()
        for chrom in chroms["chrom"]:
            chrom_length = np.array(chroms[chroms["chrom"] == chrom]["length"])[0]
            tdb = self._build_uri("genome", "sequence", f"{chrom}.tdb")
            dom = tiledb.Domain(
                tiledb.Dim(
                    name="position",
                    domain=(0, chrom_length),
                    dtype=np.int64,
                    tile=tile,
                )
            )
            attr = tiledb.Attr(
                name="nucleotide",
                dtype="ascii",
                filters=tiledb.FilterList(
                    [
                        tiledb.LZ4Filter(),
                        tiledb.ZstdFilter(level=compression),
                    ],
                    chunksize=1000,
                ),
            )
            schema = tiledb.ArraySchema(
                ctx=self.cfg.ctx,
                domain=dom,
                attrs=[attr],
                sparse=False,
                coords_filters=tiledb.FilterList(
                    [
                        tiledb.LZ4Filter(),
                        tiledb.ZstdFilter(level=compression),
                    ],
                    chunksize=1000,
                ),
            )
            tiledb.Array.create(tdb, schema)

    def _create_track_schema(self, max_bws: int, tile: int, compression: int):
        # Create /coverage/tracks.tdb
        tdb = self._build_uri("coverage", "tracks.tdb")
        dom = tiledb.Domain(
            tiledb.Dim(name="idx", domain=(0, max_bws), dtype=np.int64, tile=1),
        )
        attr1 = tiledb.Attr(name="label", dtype="ascii")
        attr2 = tiledb.Attr(name="path", dtype="ascii")
        schema = tiledb.ArraySchema(
            ctx=self.cfg.ctx, domain=dom, attrs=[attr1, attr2], sparse=False
        )
        tiledb.Array.create(tdb, schema)

        # Create every /coverage/{chrom}.tdb
        chroms = self.chroms()
        for chrom in chroms["chrom"]:
            chrom_length = np.array(chroms[chroms["chrom"] == chrom]["length"])[0]
            tdb = self._build_uri("coverage", f"{chrom}.tdb")
            dom = tiledb.Domain(
                tiledb.Dim(
                    name="position",
                    domain=(0, chrom_length),
                    dtype=np.int64,
                    tile=tile,
                )
            )
            attr = tiledb.Attr(name="placeholder", dtype="float32")
            schema = tiledb.ArraySchema(
                ctx=self.cfg.ctx,
                domain=dom,
                attrs=[attr],
                sparse=False,
                coords_filters=tiledb.FilterList(
                    [
                        tiledb.LZ4Filter(),
                        tiledb.ZstdFilter(level=compression),
                    ],
                    chunksize=1000,
                ),
            )
            tiledb.Array.create(tdb, schema)

    def _populate_track_table(self, bws: Dict[str, str]):
        n = self.tracks().shape[0]

        tdb = self._build_uri("coverage", "tracks.tdb")
        with tiledb.open(tdb, mode="w", ctx=self.cfg.ctx) as array:
            array[n : (n + len(bws))] = {
                "label": list(bws.keys()),
                "path": list(bws.values()),
            }

    def _populate_chroms_table(self, bws: Dict[str, str], threads: int):

        def _add_attribute_to_array(uri, attribute_name):
            # Check that attribute does not already exist
            has_attr = False
            with tiledb.open(uri, mode="r", ctx=self.cfg.ctx) as A:
                if A.schema.has_attr(attribute_name):
                    logger.warning(
                        f"Label {attribute_name} already exists and will be erased."
                    )
                    has_attr = True

            # Add attribute to array
            if not has_attr:
                new_attr = tiledb.Attr(attribute_name, dtype="float32")
                se = tiledb.ArraySchemaEvolution(self.cfg.ctx)
                se.add_attribute(new_attr)
                se.array_evolve(uri)

            # Check whether `placeholder` attribute still exists
            erase_placeholder = False
            with tiledb.open(uri, mode="r", ctx=self.cfg.ctx) as A:
                if A.schema.has_attr("placeholder") and A.nattr > 1:
                    erase_placeholder = True

            # If so, drop it
            if erase_placeholder:
                se = tiledb.ArraySchemaEvolution(self.cfg.ctx)
                se.drop_attribute("placeholder")
                se.array_evolve(uri)
                logger.debug(f"Attribute 'placeholder' found in {uri}. Dropping it.")

        def _process_chrom(self, chrom, chrom_length, bws):

            tdb = self._build_uri("coverage", f"{chrom}.tdb")

            # If there are already scores in the array, read them
            # THIS NEEDS TO BE UPDATED AS SOON AS
            # TILEDB ALLOWS PARTIAL ATTRIBUTE WRITING!!
            with tiledb.open(tdb, mode="r", ctx=self.cfg.ctx) as A:
                sch = A.schema
                attrs = [sch.attr(i).name for i in range(0, sch.nattr)]
                if len(attrs) == 1 and attrs[0] == "placeholder":
                    orig_scores = {}
                else:
                    orig_scores = A[0:chrom_length]

            for bwf in bws.keys():
                # Add bw label to array attributes
                _add_attribute_to_array(tdb, bwf)

                # Ingest bigwig scores for this bw and this chrom
                with pyBigWig.open(bws[bwf]) as bw:
                    arr = bw.values(chrom, 0, chrom_length, numpy=True)
                    orig_scores[bwf] = arr

            # Re-write appended scores to chrom array
            with tiledb.open(tdb, mode="w", ctx=self.cfg.ctx) as A:
                A[0:chrom_length] = orig_scores

        def _log_task_completion(future, chrom, ntasks, completed_tasks):
            if future.exception() is not None:
                logger.error(
                    f"Tracks ingestion over {chrom} failed with exception: {future.exception()}"
                )
            else:
                with lock:
                    completed_tasks[0] += 1
                logger.info(
                    f"task {completed_tasks[0]}/{ntasks} :: ingested tracks over {chrom}."
                )

        tasks = []
        chroms = self.chroms()
        for chrom in chroms["chrom"]:
            chrom_length = np.array(chroms[chroms["chrom"] == chrom]["length"])[0]
            tasks.append((chrom, chrom_length))
        ntasks = len(chroms)
        completed_tasks = [0]
        threads = min(threads, ntasks)

        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            futures = []
            for chrom, chrom_length in tasks:
                future = executor.submit(_process_chrom, self, chrom, chrom_length, bws)
                future.add_done_callback(
                    lambda f, c=chrom: _log_task_completion(
                        f, c, ntasks, completed_tasks
                    )
                )
                futures.append(future)
            concurrent.futures.wait(futures)

    def _populate_sequence_table(self, fasta: str, threads: int):

        def _process_chrom(self, chrom, chroms, fasta):
            tdb = self._build_uri("genome", "sequence", f"{chrom}.tdb")
            chrom_length = np.array(chroms[chroms["chrom"] == chrom]["length"])[0]
            with pyfaidx.Fasta(fasta) as fa:
                chrom_seq = fa.get_seq(chrom, 1, chrom_length + 1)
                chrom_seq = np.array(list(chrom_seq.seq), dtype="S1")
            with tiledb.open(tdb, mode="w", ctx=self.cfg.ctx) as A:
                A[0:chrom_length] = {"nucleotide": chrom_seq}

        def _log_task_completion(future, chrom, ntasks, completed_tasks):
            if future.exception() is not None:
                logger.error(
                    f"Fasta ingestion over {chrom} failed with exception: {future.exception()}"
                )
            else:
                with lock:
                    completed_tasks[0] += 1
                logger.info(
                    f"task {completed_tasks[0]}/{ntasks} :: ingested fasta over {chrom}."
                )

        chroms = self.chroms()
        tasks = chroms["chrom"]
        ntasks = len(tasks)
        completed_tasks = [0]
        threads = min(threads, ntasks)
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            futures = []
            for chrom in tasks:
                future = executor.submit(_process_chrom, self, chrom, chroms, fasta)
                future.add_done_callback(
                    lambda f, c=chrom: _log_task_completion(
                        f, c, ntasks, completed_tasks
                    )
                )
                futures.append(future)
            concurrent.futures.wait(futures)

    def chroms(self) -> pd.DataFrame:
        """Extract chromosome table from a `.momics` repository.

        Returns:
            pd.DataFrame: A data frame listing one chromosome per row
        """
        try:
            chroms = self._get_table(self._build_uri("genome", "chroms.tdb"))
        except FileExistsError:
            chroms = pd.DataFrame(columns=["chrom_index", "chrom", "length"])
        return chroms

    def seq(self) -> pd.DataFrame:
        """Extract sequence table from a `.momics` repository.

        Returns:
            pd.DataFrame: A data frame listing one chromosome per row, with first/last 10 nts.
        """
        if self.chroms().empty:
            raise OSError("`chroms` table has not been filled out yet.")

        try:
            tdb = self._build_uri(
                "genome", "sequence", f"{self.chroms()['chrom'][0]}.tdb"
            )
            _ = self._get_table(tdb)
            pass
        except FileExistsError:
            raise OSError("`seq` table has not been filled out yet.")

        chroms = self.chroms()
        chroms["seq"] = pd.Series()
        for chrom in chroms["chrom"]:
            tdb = self._build_uri("genome", "sequence", f"{chrom}.tdb")
            chrom_len = chroms[chroms["chrom"] == chrom]["length"].iloc[0]
            with tiledb.open(tdb, "r", ctx=self.cfg.ctx) as A:
                start_nt = "".join(A.df[0:9]["nucleotide"])
                end_nt = "".join(A.df[(chrom_len - 10) : (chrom_len - 1)]["nucleotide"])
            chroms.loc[chroms["chrom"] == chrom, "seq"] = start_nt + "..." + end_nt

        return chroms

    def tracks(self) -> pd.DataFrame:
        """Extract table of ingested bigwigs.

        Returns:
            pd.DataFrame: A data frame listing one ingested bigwig file per row
        """
        try:
            tracks = self._get_table(self._build_uri("coverage", "tracks.tdb"))
            tracks = tracks[tracks["label"] != "\x00"]
        except FileExistsError:
            tracks = pd.DataFrame(columns=["idx", "label", "path"])
        return tracks

    def bins(self, width, step, cut_last_bin_out=False):
        """Generate a DataFrame of tiled genomic bins

        Args:
            width (_type_): The width of each bin.
            step (_type_): The step size for tiling.
            cut_last_bin_out (bool, optional): Remove the last bin of each chromosome. Defaults to False.

        Returns:
            _type_: pd.DataFrame: DataFrame with columns "chrom", 'start', 'end'.
        """
        bins = []
        chroms = self.chroms().set_index("chrom")["length"].to_dict()

        for chrom, length in chroms.items():
            start = 0
            while start < length:
                end = min(start + width, length)
                bins.append({"chrom": chrom, "start": (start + 1), "end": end})
                start += step

        df = pd.DataFrame(bins)
        if cut_last_bin_out:
            df = df[(df["end"] - df["start"]) == width - 1]

        return df

    def add_chroms(self, chr_lengths: dict, genome_version: str = "") -> "Momics":
        """Add chromosomes (and genome) information the `.momics` repository.

        Args:
            chr_lengths (dict): Chromosome lengths
            genome_version (str, optional): Genome version (default: ""). Defaults to "".

        Returns:
            Momics: An updated Momics object
        """
        if not self.chroms().empty:
            raise ValueError("`chroms` table has already been filled out.")

        tdb = self._build_uri("genome", "chroms.tdb")
        dom_genome = tiledb.Domain(
            tiledb.Dim(
                name="chrom_index",
                domain=(0, len(chr_lengths) - 1),
                dtype=np.int32,
                tile=len(chr_lengths),
            )
        )
        attr_chr = tiledb.Attr(name="chrom", dtype="ascii", var=True)
        attr_length = tiledb.Attr(name="length", dtype=np.int64)
        schema = tiledb.ArraySchema(
            ctx=self.cfg.ctx,
            domain=dom_genome,
            attrs=[attr_chr, attr_length],
            sparse=False,
        )
        tiledb.Array.create(tdb, schema)

        # Populate `chrom` array
        chr = list(chr_lengths.keys())
        length = list(chr_lengths.values())
        with tiledb.open(tdb, "w", ctx=self.cfg.ctx) as A:
            A[0 : len(chr)] = {"chrom": np.array(chr, dtype="S"), "length": length}
            A.meta["genome_assembly_version"] = genome_version
            A.meta["timestamp"] = datetime.now().isoformat()

    def add_sequence(
        self, fasta: Path, threads: int = 1, tile: int = 10000, compression: int = 3
    ) -> "Momics":
        """Ingest a fasta file into a Momics repository

        Args:
            fasta (str): Path to a Fasta file containing the genome reference sequence.
            threads (int, optional): Threads to parallelize I/O. Defaults to 1.
            tile (int, optional): Tile size for TileDB. Defaults to 10000.
            compression (int, optional): Compression level for TileDB. Defaults to 3.

        Returns:
            Momics: The updated Momics object
        """
        # Abort if `chroms` have not been filled
        chroms = self.chroms()
        if chroms.empty:
            raise ValueError("Please fill out `chroms` table first.")

        # Abort if sequence table already exists
        tdb = self._build_uri("genome", "sequence", f"{chroms['chrom'][0]}.tdb")
        if self.cfg.vfs.is_dir(tdb):
            raise tiledb.cc.TileDBError(f"Error: TileDB '{tdb}' already exists.")

        # Abort if chr lengths in provided fasta do not match those in `chroms`
        utils._check_fasta_lengths(fasta, chroms)

        # Create sequence tables schema
        self._create_sequence_schema(tile, compression)

        # Populate each `/genome/sequence/{chrom}.tdb`
        self._populate_sequence_table(fasta, threads)

    def add_tracks(
        self,
        bws: dict,
        threads: int = 1,
        max_bws: int = 9999,
        tile: int = 10000,
        compression: int = 3,
    ) -> "Momics":
        """Ingest bigwig coverage tracks to the `.momics` repository.

        Args:
            bws (dict): Dictionary of bigwig files
            threads (int, optional): Threads to parallelize I/O. Defaults to 1.
            max_bws (int, optional): Maximum number of bigwig files. Defaults to 9999.
            tile (int, optional): Tile size. Defaults to 10000.
            compression (int, optional): Compression level. Defaults to 3.

        Returns:
            Momics: The updated Momics object
        """
        start0 = time.time()

        # Abort if `chroms` have not been filled
        if self.chroms().empty:
            raise ValueError("Please fill out `chroms` table first.")

        # Abort if chr lengths in provided bw do not match those in `chroms`
        utils._check_chr_lengths(bws, self.chroms())

        # Abort if bw labels already exist
        utils._check_track_names(bws, self.tracks())

        # If `path/coverage/tracks.tdb` (and `{chroms.tdb}`) do not exist, create it
        if self.tracks().empty:
            self._create_track_schema(max_bws, tile, compression)

        # Populate each `path/coverage/{chrom}.tdb`
        self._populate_chroms_table(bws, threads)

        # Populate `path/coverage/tracks.tdb`
        self._populate_track_table(bws)

        logger.info(f"{len(bws)} tracks ingested in {round(time.time() - start0,4)}s.")

    def add_track(
        self,
        coverage: dict,
        track: str,
        threads: int = 1,
    ) -> "Momics":
        """
        Ingest a coverage track provided as a dictionary to a `.momics` repository.
        This method is useful when you have already computed the coverage track and
        have it in memory.

        Args:
            coverage (dict): Dictionary of coverage tracks. The keys are chromosome names and the values are numpy arrays.
            track (str): Label to store the track under.
            threads (int, optional): Threads to parallelize I/O. Defaults to 1.

        Returns:
            Momics: The updated Momics object
        """
        chroms = self.chroms()
        tracks = self.tracks()

        # Abort if `chroms` have not been filled
        if chroms.empty:
            raise ValueError("Please fill out `chroms` table first.")
        if tracks.empty:
            raise ValueError("Please fill out `tracks` table first.")

        # Abort if chr lengths in provided bw do not match those in `chroms`
        reference_lengths = dict(zip(chroms["chrom"], chroms["length"]))
        lengths = dict(zip(chroms["chrom"], [len(v) for k, v in coverage.items()]))
        if lengths != reference_lengths:
            raise Exception(
                f"`{track}` coverage track does not chromomosome lengths matching those of the momics repository."
            )

        # Abort if bw labels already exist
        if track in set(tracks["label"]):
            raise ValueError(
                f"Provided label '{track}' already present in `tracks` table"
            )

        # Save the coverage dict as a temporary bigwig file
        # and ingest it using `add_tracks`
        tmp_bw = tempfile.NamedTemporaryFile(delete=False)
        utils._dict_to_bigwig(coverage, tmp_bw.name)
        self.add_tracks({track: tmp_bw.name}, threads=threads)
        os.remove(tmp_bw.name)

    def remove_track(self, track: str) -> "Momics":
        """Remove a track from a `.momics` repository.

        Args:
            track (str): Which track to remove

        Returns:
            Momics: An updated Momics object
        """
        # Abort if `track` is not listed
        tracks = self.tracks()
        chroms = self.chroms()
        utils._check_track_name(track, tracks)

        # Remove entry from each `path/coverage/{chrom}.tdb`
        # and from `path/coverage/tracks.tdb`
        for chrom in chroms["chrom"]:
            tdb = self._build_uri("coverage", f"{chrom}.tdb")
            ctx = self.cfg.ctx
            se = tiledb.ArraySchemaEvolution(ctx)
            se.drop_attribute(track)
            se.array_evolve(tdb)

        tdb = self._build_uri("coverage", "tracks.tdb")
        idx = tracks["idx"][tracks["label"] == track].values[0]
        with tiledb.open(tdb, mode="w", ctx=self.cfg.ctx) as A:
            A[idx] = {"label": None, "path": None}

    def remove(self) -> bool:
        """Remove a `.momics` repository."""
        host = self._is_cloud_hosted()
        vfs = self.cfg.vfs

        ## Remove local repo
        if not host:
            vfs.remove_dir(self.path)
            logger.info(f"Purged {self.path}")

        ## Remove S3 and GCS-hosted repo
        if host in ["s3", "gcs"]:
            vfs.remove_dir(self.path)
            logger.info(f"Purged {self.path}")

        if host == "azure":

            def remove_directory_until_success(
                vfs, dir_uri, max_retries=10, retry_delay=2
            ):
                attempts = 0
                while attempts < max_retries:
                    try:
                        vfs.remove_dir(dir_uri)
                        logger.info(f"Purged {dir_uri}")
                        break
                    except tiledb.TileDBError as e:
                        attempts += 1
                        if attempts < max_retries:
                            time.sleep(retry_delay)
                        else:
                            raise e

            remove_directory_until_success(vfs, self.path)

        return True

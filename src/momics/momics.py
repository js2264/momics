from datetime import datetime
from time import sleep
from typing import Dict, Optional
from pathlib import Path
import concurrent.futures
from .logging import logger
import threading

import numpy as np
import pandas as pd
import pyBigWig
import pyfaidx
import tiledb

from . import utils
from .config import MomicsConfig


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
        create=True,
        config: Optional[MomicsConfig] = None,
    ):
        """Initialize the Momics class.

        Args:
            path (str): Path to a `.momics` repository.
            create (bool, optional): If not found, should the repository be initiated? Defaults to True.
        """
        self.path = path
        if config is None:
            config = MomicsConfig()
        self.cfg = config

        ## Check if folder exists:
        if self.cfg.vfs.is_dir(self.path):
            ## If it exists and `create` = True, raise Error
            if create:
                raise OSError(f"{path} repository already exists.")
        ## If does not exist:
        else:
            ## Check if we want to create it
            if create:
                self._create_repository()
            else:
                raise OSError("Momics repository not found.")

    def _is_cloud_hosted(self):
        if self.path.startswith(("s3://", "gcs://", "azure://")):
            return True

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
            a = A.df[:]

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
                    domain=(0, chrom_length - 1),
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
            tiledb.DenseArray.create(tdb, schema)

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
                    domain=(0, chrom_length - 1),
                    dtype=np.int64,
                    tile=tile,
                ),
                tiledb.Dim(name="idx", domain=(0, max_bws), dtype=np.int64, tile=1),
            )
            attr = tiledb.Attr(
                name="scores",
                dtype=np.float32,
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
                sparse=True,
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
        try:
            n = self.tracks().shape[0]
        except tiledb.cc.TileDBError:
            n = 0

        tdb = self._build_uri("coverage", "tracks.tdb")
        with tiledb.DenseArray(tdb, mode="w", ctx=self.cfg.ctx) as array:
            array[n : (n + len(bws))] = {
                "label": list(bws.keys()),
                "path": list(bws.values()),
            }

    def _populate_chroms_table(self, bws: Dict[str, str], threads: int):
        try:
            n = self.tracks().shape[0]
        except tiledb.cc.TileDBError:
            n = 0

        def _process_chrom_bwf(self, chrom, chrom_length, idx, bwf, bws, n):
            sleep(1)
            tdb = self._build_uri("coverage", f"{chrom}.tdb")
            with pyBigWig.open(bws[bwf]) as bw:
                arr = np.array(bw.values(chrom, 0, chrom_length), dtype=np.float32)
                with tiledb.open(tdb, mode="w", ctx=self.cfg.ctx) as A:
                    coord1 = np.arange(0, chrom_length)
                    coord2 = np.repeat(idx + n, len(coord1))
                    A[coord1, coord2] = {"scores": arr}

        def _log_task_completion(future, bwf, chrom, ntasks, completed_tasks):
            if future.exception() is not None:
                logger.error(
                    f"Bigwig {bwf} ingestion over {chrom} failed with exception: {future.exception()}"
                )
            else:
                with lock:
                    completed_tasks[0] += 1
                logger.info(
                    f"task {completed_tasks[0]}/{ntasks} :: ingested bigwig {bwf} over {chrom}."
                )

        tasks = []
        chroms = self.chroms()
        for chrom in chroms["chrom"]:
            chrom_length = np.array(chroms[chroms["chrom"] == chrom]["length"])[0]
            for idx, bwf in enumerate(bws):
                tasks.append((chrom, chrom_length, idx, bwf))
        ntasks = len(tasks)
        completed_tasks = [0]
        threads = min(threads, ntasks)

        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            futures = []
            for chrom, chrom_length, idx, bwf in tasks:
                future = executor.submit(
                    _process_chrom_bwf, self, chrom, chrom_length, idx, bwf, bws, n
                )
                future.add_done_callback(
                    lambda f, b=bwf, c=chrom: _log_task_completion(
                        f, b, c, ntasks, completed_tasks
                    )
                )
                futures.append(future)
            concurrent.futures.wait(futures)

    def _populate_sequence_table(self, fasta: str, threads: int):

        def _process_chrom(self, chrom, chroms, fasta):
            sleep(1)
            tdb = self._build_uri("genome", "sequence", f"{chrom}.tdb")
            chrom_length = np.array(chroms[chroms["chrom"] == chrom]["length"])[0]
            with pyfaidx.Fasta(fasta) as fa:
                chrom_seq = fa.get_seq(chrom, 1, chrom_length)
            with tiledb.DenseArray(tdb, mode="w", ctx=self.cfg.ctx) as A:
                A[:] = {"nucleotide": np.array(list(chrom_seq.seq), dtype="S1")}

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
            sparse=True,
        )
        tiledb.Array.create(tdb, schema)

        # Populate `chrom` array
        chr = list(chr_lengths.keys())
        length = list(chr_lengths.values())
        with tiledb.open(tdb, "w", ctx=self.cfg.ctx) as array:
            indices = np.arange(len(chr))
            array[indices] = {"chrom": np.array(chr, dtype="S"), "length": length}
            array.meta["genome_assembly_version"] = genome_version
            array.meta["timestamp"] = datetime.now().isoformat()

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

    def remove_track(self, track: str) -> "Momics":
        """Remove a track from a `.momics` repository.

        Args:
            track (str): Which track to remove

        Returns:
            Momics: An updated Momics object
        """
        # Abort if `track` is not listed
        utils._check_track_name(track, self.tracks())

        # Remove entry from each `path/coverage/{chrom}.tdb`
        # and from `path/coverage/tracks.tdb`
        idx = self.tracks()["idx"][self.tracks()["label"] == track].values[0]
        qc = f"idx == {idx}"
        for chrom in self.chroms()["chrom"]:
            tdb = self._build_uri("coverage", f"{chrom}.tdb")
            with tiledb.open(tdb, mode="d", ctx=self.cfg.ctx) as A:
                A.query(cond=qc).submit()

        tdb = self._build_uri("coverage", "tracks.tdb")
        with tiledb.open(tdb, mode="w", ctx=self.cfg.ctx) as A:
            A[idx] = {"label": None, "path": None}

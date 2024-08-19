import os
from datetime import datetime
from typing import Dict, Optional

import numpy as np
import pandas as pd
import pyBigWig
import tiledb

from . import utils


class Momics:
    """
    A class to manipulate `.momics` repositories.

    Attributes
    ----------
    path : str
        Path to a `.momics` repository.
    """

    def __init__(self, path: str, create=True):
        """
        Initialize the Momics class.

        Parameters
        ----------
        path : str
            Path to a `.momics` repository.

        create : bool
            If not found, should the repository be initiated?
        """

        self.path = path

        if not os.path.exists(path):
            if create:
                self._create_repository()
            else:
                raise OSError("Momics repository not found.")
        else:
            if create:
                raise OSError(f"{path} repository already exist found.")

    def _create_repository(self):
        genome_path = os.path.join(self.path, "genome")
        seq_path = os.path.join(self.path, "genome", "sequence")
        coverage_path = os.path.join(self.path, "coverage")
        features_path = os.path.join(self.path, "features")
        tiledb.group_create(self.path)
        tiledb.group_create(genome_path)
        tiledb.group_create(coverage_path)
        tiledb.group_create(features_path)
        tiledb.group_create(seq_path)

    def _get_table(self, tdb: str) -> Optional[pd.DataFrame]:
        path = os.path.join(self.path, tdb)
        if not os.path.exists(path):
            return None

        with tiledb.open(path, "r") as A:
            a = A.df[:]

        return a

    def _create_track_schema(self, max_bws: int, tile: int, compression: int):
        # Create path/coverage/tracks.tdb
        tdb = os.path.join(self.path, "coverage", "tracks.tdb")
        dom = tiledb.Domain(
            tiledb.Dim(name="idx", domain=(0, max_bws), dtype=np.int64, tile=1),
        )
        attr1 = tiledb.Attr(name="label", dtype="ascii")
        attr2 = tiledb.Attr(name="path", dtype="ascii")
        schema = tiledb.ArraySchema(domain=dom, attrs=[attr1, attr2], sparse=False)
        tiledb.Array.create(tdb, schema)
        chroms = self.chroms()

        # Create every path/coverage/{chrom}.tdb
        for chrom in chroms["chr"]:
            chrom_length = np.array(chroms[chroms["chr"] == chrom]["length"])[0]
            tdb = os.path.join(self.path, "coverage", f"{chrom}.tdb")
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

        tdb = os.path.join(self.path, "coverage", "tracks.tdb")
        with tiledb.DenseArray(tdb, mode="w") as array:
            array[n : (n + len(bws))] = {
                "label": list(bws.keys()),
                "path": list(bws.values()),
            }

    def _populate_chroms_table(self, bws: Dict[str, str]):
        try:
            n = self.tracks().shape[0]
        except tiledb.cc.TileDBError:
            n = 0

        chroms = self.chroms()
        for chrom in chroms["chr"]:
            chrom_length = np.array(chroms[chroms["chr"] == chrom]["length"])[0]
            tdb = os.path.join(self.path, "coverage", f"{chrom}.tdb")
            for idx, bwf in enumerate(bws):
                with pyBigWig.open(bws[bwf]) as bw:
                    arr = np.array(bw.values(chrom, 0, chrom_length), dtype=np.float32)
                    with tiledb.open(tdb, mode="w") as A:
                        coord1 = np.arange(0, chrom_length)
                        coord2 = np.repeat(idx + n, len(coord1))
                        A[coord1, coord2] = {"scores": arr}

    def chroms(self) -> pd.DataFrame:
        """
        Extract chromosome table from a `.momics` repository.

        Returns
        -------
        pd.DataFrame
            A data frame listing one chromosome per row
        """
        chroms = self._get_table(os.path.join("genome", "chroms.tdb"))
        return (
            chroms
            if chroms is not None
            else pd.DataFrame(columns=["chrom_index", "chr", "length"])
        )

    def tracks(self) -> pd.DataFrame:
        """
        Extract table of ingested bigwigs.

        Returns
        -------
        pandas.DataFrame
            A data frame listing one ingested bigwig file per row
        """
        tracks = self._get_table(os.path.join("coverage", "tracks.tdb"))
        return (
            tracks
            if tracks is not None
            else pd.DataFrame(columns=["idx", "label", "path"])
        )

    def add_tracks(
        self, bws: dict, max_bws: int = 9999, tile: int = 10000, compression: int = 3
    ):
        """
        Ingest bigwig coverage tracks to the `.momics` repository.

        Parameters
        ----------
        bws : dict
            Dictionary of bigwig files
        max_bws : int
            Maximum number of bigwig files.
        tile : int
            Tile size.
        compression : int
            Compression level.
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
        self._populate_chroms_table(bws)

        # Populate `path/coverage/tracks.tdb`
        self._populate_track_table(bws)

    def add_chroms(self, chr_lengths: dict, genome_version: str = ""):
        """
        Add chromosomes (and genome) information the `.momics` repository.

        Parameters
        ----------
        chr_lengths : dict
            Chromosome lengths
        genome_version : str
            Genome version (default: "")
        """

        if not self.chroms().empty:
            raise ValueError("`chroms` table has already been filled out.")

        tdb = os.path.join(self.path, "genome", "chroms.tdb")
        dom_genome = tiledb.Domain(
            tiledb.Dim(
                name="chrom_index",
                domain=(0, len(chr_lengths) - 1),
                dtype=np.int32,
                tile=len(chr_lengths),
            )
        )
        attr_chr = tiledb.Attr(name="chr", dtype="ascii", var=True)
        attr_length = tiledb.Attr(name="length", dtype=np.int64)
        schema = tiledb.ArraySchema(
            domain=dom_genome, attrs=[attr_chr, attr_length], sparse=True
        )
        tiledb.Array.create(tdb, schema)

        # Populate `chrom` array
        chr = list(chr_lengths.keys())
        length = list(chr_lengths.values())
        with tiledb.open(tdb, "w") as array:
            indices = np.arange(len(chr))
            array[indices] = {"chr": np.array(chr, dtype="S"), "length": length}
            array.meta["genome_assembly_version"] = genome_version
            array.meta["timestamp"] = datetime.now().isoformat()

    def query(
        self,
        query: str,
    ):
        """
        Query bigwig coverage tracks from a `.momics` repository.

        Parameters
        ----------
        chr_lengths : str
            UCSC-style chromosome interval (e.g. "II:12001-15000")
        """
        tr = self.tracks().drop("path", axis=1)
        chrom, range_part = query.split(":")
        start = range_part.split("-")[0]
        end = range_part.split("-")[1]
        with tiledb.open(f"test.momics/coverage/{chrom}.tdb", "r") as array:
            data = array.df[int(start) : int(end), :]

        data["idx"] = data["idx"] - 1
        mdata = pd.merge(tr, data, on="idx").drop("idx", axis=1)
        return mdata

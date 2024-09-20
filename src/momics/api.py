import os
from datetime import datetime
from typing import Dict, Optional

import numpy as np
import pandas as pd
import pyBigWig
import pyfaidx
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

    def _create_sequence_schema(self, tile: int, compression: int):
        # Create every path/genome/sequence/{chrom}.tdb
        chroms = self.chroms()
        for chrom in chroms["chr"]:
            chrom_length = np.array(chroms[chroms["chr"] == chrom]["length"])[0]
            tdb = os.path.join(self.path, "genome", "sequence", f"{chrom}.tdb")
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
        # Create path/coverage/tracks.tdb
        tdb = os.path.join(self.path, "coverage", "tracks.tdb")
        dom = tiledb.Domain(
            tiledb.Dim(name="idx", domain=(0, max_bws), dtype=np.int64, tile=1),
        )
        attr1 = tiledb.Attr(name="label", dtype="ascii")
        attr2 = tiledb.Attr(name="path", dtype="ascii")
        schema = tiledb.ArraySchema(domain=dom, attrs=[attr1, attr2], sparse=False)
        tiledb.Array.create(tdb, schema)

        # Create every path/coverage/{chrom}.tdb
        chroms = self.chroms()
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

    def _populate_sequence_table(self, fasta: str):
        chroms = self.chroms()
        for chrom in chroms["chr"]:
            tdb = os.path.join(self.path, "genome", "sequence", f"{chrom}.tdb")
            chrom_length = np.array(chroms[chroms["chr"] == chrom]["length"])[0]
            with pyfaidx.Fasta(fasta) as fa:
                chrom_seq = fa.get_seq(chrom, 1, chrom_length)
            with tiledb.DenseArray(tdb, mode="w") as A:
                A[:] = {"nucleotide": np.array(list(chrom_seq.seq), dtype="S1")}

    def _purge_track(self, track: str):
        idx = self.tracks()["idx"][self.tracks()["label"] == track].values[0]
        qc = f"idx == {idx}"
        for chrom in self.chroms()["chr"]:
            tdb = os.path.join(self.path, "coverage", f"{chrom}.tdb")
            with tiledb.open(tdb, mode="d") as A:
                A.query(cond=qc).submit()

        tdb = os.path.join(self.path, "coverage", "tracks.tdb")
        with tiledb.open(tdb, mode="w") as A:
            A[idx] = {"label": None, "path": None}

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

    def seq(self) -> pd.DataFrame:
        """
        Extract sequence table from a `.momics` repository.

        Returns
        -------
        pd.DataFrame
            A data frame listing one chromosome per row, with first/last 10 nts.
        """
        try:
            self.query_sequence(f"{self.chroms()['chr'][0]}:1-2")
        except tiledb.cc.TileDBError:
            raise ValueError("Genomic sequence not added yet to the repository.")

        chroms = self.chroms()
        chroms["seq"] = pd.Series()
        for chrom in chroms["chr"]:
            chrom_len = chroms[chroms["chr"] == chrom]["length"].iloc[0]
            start_nt = self.query_sequence(f"{chrom}:1-10")
            end_nt = self.query_sequence(f"{chrom}:{chrom_len-10}-{chrom_len}")
            chroms.loc[chroms["chr"] == chrom, "seq"] = start_nt + "..." + end_nt

        return chroms

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
            # tracks[tracks["label"] != "None"]
            tracks
            if tracks is not None
            else pd.DataFrame(columns=["idx", "label", "path"])
        )

    def bins(self, width, step, cut_last_bin_out=True):
        """
        Generate a DataFrame of tiled genomic bins.

        Parameters:
        width (int): The width of each bin.
        step (int): The step size for tiling.

        Returns:
        pd.DataFrame: DataFrame with columns 'chr', 'start', 'end'.
        """

        bins = []
        chroms = self.chroms().set_index("chr")["length"].to_dict()

        for chrom, length in chroms.items():
            start = 0
            while start < length:
                end = min(start + width, length)
                bins.append({"chr": chrom, "start": (start + 1), "end": end})
                start += step

        df = pd.DataFrame(bins)
        if cut_last_bin_out:
            df = df[(df["end"] - df["start"]) == width - 1]

        return df

    def add_sequence(self, fasta: str, tile: int = 10000, compression: int = 3):
        """
        Ingest multi-sequence fasta file to the `.momics` repository.

        Parameters
        ----------
        fasta : str
            Path to fasta file
        tile : int
            Tile size.
        compression : int
            Compression level.
        """

        # Abort if `chroms` have not been filled
        if self.chroms().empty:
            raise ValueError("Please fill out `chroms` table first.")

        # Abort if sequence table already exists
        try:
            self.query_sequence(f"{self.chroms()['chr'][0]}:1-2")
            raise ValueError("Sequence already added to the repository.")
        except tiledb.cc.TileDBError:
            pass

        # Abort if chr lengths in provided fasta do not match those in `chroms`
        utils._check_fasta_lengths(fasta, self.chroms())

        # Create sequence tables schema
        self._create_sequence_schema(tile, compression)

        # Populate each `path/genome/sequence/{chrom}.tdb`
        self._populate_sequence_table(fasta)

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

    def remove_track(self, track: str):
        """
        Remove a track from a `.momics` repository.

        Parameters
        ----------
        track : str
            Which track to remove
        """

        # Abort if `track` is not listed
        utils._check_track_name(track, self.tracks())

        # Remove entry from each `path/coverage/{chrom}.tdb` and from `path/coverage/tracks.tdb`
        self._purge_track(track)

    def export_track(self, track: str, prefix: str):
        """
        Export a track from a `.momics` repository as a .bw file.

        Parameters
        ----------
        track : str
            Which track to remove
        prefix : str
            Prefix of the output bigwig file
        """

        # Abort if `track` is not listed
        utils._check_track_name(track, self.tracks())

        # Init output file
        output = prefix + ".bw"
        bw = pyBigWig.open(output, "w")
        chrom_sizes = self.chroms()[["chr", "length"]].apply(tuple, axis=1).tolist()
        bw.addHeader(chrom_sizes)
        for chrom, _ in chrom_sizes:
            print(chrom)
            q = self.query_tracks(chrom)
            q = q[q["label"] == track].dropna(subset=["scores"])
            # starts = q["position"].to_numpy(dtype=np.int64)
            # ends = starts + 1
            # values = q["scores"].to_numpy(dtype=np.int64)
            # chroms = np.array([chrom] * 10)

            starts = q["position"].to_numpy(dtype=np.int64)
            ends = starts + 1
            values0 = q["scores"].to_numpy(dtype=np.float32)
            chroms = np.array([chrom] * len(values0))
            bw.addEntries(chroms, starts=starts, ends=ends, values=values0)
        bw.close()

    def query_sequence(
        self,
        query: str,
    ):
        """
        Query chromosome sequence from a `.momics` repository.

        Parameters
        ----------
        chr_lengths : str
            UCSC-style chromosome interval (e.g. "II:12001-15000")
        """
        if ":" in query:
            chrom, range_part = query.split(":")
            utils._check_chr_name(chrom, self.chroms())
            start = int(range_part.split("-")[0]) - 1
            end = int(range_part.split("-")[1]) - 1
            with tiledb.open(f"{self.path}/genome/sequence/{chrom}.tdb", "r") as A:
                seq = A.df[start:end]["nucleotide"]
        else:
            chrom = query
            utils._check_chr_name(chrom, self.chroms())
            with tiledb.open(f"{self.path}/genome/sequence/{chrom}.tdb", "r") as A:
                seq = A.df[:]["nucleotide"]

        return "".join(seq)

import pandas as pd
from momics import utils
import os
import tiledb
import numpy as np
import pyBigWig
from datetime import datetime


class Momics:
    """
    A class to manipulate `.momics` repositories.

    Attributes
    ----------
    path : str
        Path to a `.momics` repository.
    """

    def _get_table(self, tdb: str):
        path = os.path.join(self.path, tdb)
        if os.path.exists(path) == False:
            return None

        with tiledb.open(path, "r") as A:
            a = A.df[:]

        return a

    def chroms(self):
        """
        A method to extract chromosome table from a `.momics` repository.

        Returns
        -------
        pandas.DataFrame
            A data frame listing one chromosome per row
        """
        tdb = os.path.join("genome", "chroms.tdb")
        chroms = self._get_table(tdb)
        if chroms is None:
            return pd.DataFrame(columns=["chrom_index", "chr", "length"])

        return chroms

    def tracks(self):
        """
        Extract table of ingested bigwigs

        Returns
        -------
        pandas.DataFrame
            A data frame listing one ingested bigwig file per row
        """
        tdb = os.path.join("coverage", "tracks.tdb")
        tracks = self._get_table(tdb)
        if tracks is None:
            return pd.DataFrame(columns=["idx", "label", "path"])

        return tracks

    def add_tracks(self, bws: dict):
        """
        A method to ingest big wig coverage tracks to the `.momics` repository.

        Parameters
        ----------
        bws : dict
            Dictionary of bigwig files
        """

        # Abort if `chroms` have not been filled
        if self.chroms().empty == True:
            raise ValueError("Please fill out `chroms` table first.")

        # Abort if chr lengths in provided bw do not match those in `chroms`
        utils._check_chr_lengths(bws, self.chroms())

        # Abort if bw labels already exist
        utils._check_track_names(bws, self.tracks())

        # If `path/coverage/tracks.tdb` (and `{chroms.tdb}`) do not exist, create it
        tdb = os.path.join(self.path, "coverage", "tracks.tdb")
        if self.tracks().empty == True:
            dom = tiledb.Domain(
                tiledb.Dim(name="idx", domain=(0, 9999), dtype=np.int64),
            )
            attr1 = tiledb.Attr(name="label", dtype="ascii")
            attr2 = tiledb.Attr(name="path", dtype="ascii")
            schema = tiledb.ArraySchema(domain=dom, attrs=[attr1, attr2], sparse=False)
            tiledb.Array.create(tdb, schema)
            chroms = self.chroms()
            for chrom in chroms["chr"]:
                chrom_length = np.array(chroms[chroms["chr"] == chrom]["length"])[0]
                tdb = os.path.join(self.path, "coverage", f"{chrom}.tdb")
                dom = tiledb.Domain(
                    tiledb.Dim(
                        name="position",
                        domain=(0, chrom_length - 1),
                        dtype=np.int64,
                    ),
                    tiledb.Dim(
                        name="idx",
                        domain=(0, 999),
                        dtype=np.int64,
                    ),
                )
                attr = tiledb.Attr(
                    name="scores",
                    dtype=np.float32,
                    filters=tiledb.FilterList(
                        [tiledb.ZstdFilter(level=-3)], chunksize=10000
                    ),
                )
                schema = tiledb.ArraySchema(
                    domain=dom,
                    attrs=[attr],
                    sparse=True,
                    coords_filters=tiledb.FilterList(
                        [tiledb.ZstdFilter(level=-3)], chunksize=10000
                    ),
                )
                tiledb.Array.create(tdb, schema)
            n = 0
        else:
            n = self.tracks().shape[0]

        # Populate `path/coverage/tracks.tdb`
        tdb = os.path.join(self.path, "coverage", "tracks.tdb")
        with tiledb.DenseArray(tdb, mode="w") as array:
            array[n : (n + len(bws))] = {
                "label": list(bws.keys()),
                "path": list(bws.values()),
            }

        # Populate each `path/coverage/{chrom}.tdb`
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

    def add_chroms(self, chr_lengths: dict, genome_version: str = ""):
        """
        A method to add chromosomes (and genome) information the `.momics` repository.

        Parameters
        ----------
        chr_lengths : dict
            Chromosome lengths
        genome_version : str
            Genome version (default: "")
        """

        if self.chroms().empty == False:
            raise ValueError("`chroms` table has already been filled out.")

        tdb = os.path.join(self.path, "genome", "chroms.tdb")
        dom_genome = tiledb.Domain(
            tiledb.Dim(
                name="chrom_index", domain=(0, len(chr_lengths) - 1), dtype=np.int32
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

    def __init__(self, path: str):
        """
        Initialize the Momics class.

        Parameters
        ----------
        path : str
            Path to a `.momics` repository.
        """

        self.path = path
        genome_path = os.path.join(path, "genome")
        seq_path = os.path.join(path, "genome", "sequence")
        coverage_path = os.path.join(path, "coverage")
        features_path = os.path.join(path, "features")

        if os.path.exists(path) == False:
            tiledb.group_create(path)
            tiledb.group_create(genome_path)
            tiledb.group_create(coverage_path)
            tiledb.group_create(features_path)
            tiledb.group_create(seq_path)

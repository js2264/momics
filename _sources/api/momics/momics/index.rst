momics.momics
=============

.. py:module:: momics.momics


Attributes
----------

.. autoapisummary::

   momics.momics.BULK_OVERWRITE
   momics.momics.TILEDB_CHUNKSIZE
   momics.momics.TILEDB_COMPRESSION
   momics.momics.TILEDB_COV_FILTERS
   momics.momics.TILEDB_POSITION_FILTERS
   momics.momics.TILEDB_SEQ_FILTERS
   momics.momics.lock


Classes
-------

.. autoapisummary::

   momics.momics.Momics


Module Contents
---------------

.. py:class:: Momics(path, config = None)

   A class to manipulate `.momics` repositories.

   `.momics` repositories are a TileDB-backed storage system for genomics data.
   They are structured as follows:

   - `./genome/chroms.tdb` - table for ingested chromosomes;
   - `./coverage/tracks.tdb` - table for ingested bigwig tracks;
   - `./annotations/features.tdb` - table for ingested feature sets.

   In each subdirectory, there is also one `.tdb` file per chromosome, which
   stores the following data:

   - In `./genome/{X}.tdb`: the reference sequence of the chromosome;
   - In `./coverage/{X}.tdb`: the coverage scores of the chromosome;
   - In `./annotations/{X}.tdb`: the genomic features of the chromosome.

   .. attribute:: path

      Path to a `.momics` repository.

      :type: str

   .. attribute:: cfg

      Configuration object.

      :type: MomicsConfig


   .. py:method:: bins(width, stride = None, cut_last_bin_out=False)

      Generate a PyRanges of tiled genomic bins

      :param width: The width of each bin.
      :type width: _type_
      :param stride: The stride size for tiling.
      :type stride: _type_
      :param cut_last_bin_out: Remove the last bin of each                 chromosome. Defaults to False.
      :type cut_last_bin_out: bool, optional
      :param Remember that PyRanges are 0-based and end-exclusive.:

      :returns: pr.PyRanges: a PyRanges object of tiled genomic bins.
      :rtype: _type_



   .. py:method:: chroms(as_dict = False)

      Extract chromosome table from a `.momics` repository.

      :param as_dict: Whether to return chromosomes as a dictionary. Defaults to False.
      :type as_dict: bool, optional

      :returns:     - If as_dict=False: A data frame listing one chromosome per row
                    - If as_dict=True: A dictionary with chromosome names as keys and lengths as values
      :rtype: Union[pd.DataFrame, Dict[str, int]]



   .. py:method:: consolidate(vacuum = True)

      Consolidates the fragments of all arrays in the repository.

      :param vacuum: Vacuum the consolidated array. Defaults to True.
      :type vacuum: bool, optional



   .. py:method:: export_features(features, output)

      Export a features set from a `.momics` repository as a `.bed` file.

      :param features: Which features to remove
      :type features: str
      :param output: Prefix of the output BED file
      :type output: Path

      :returns: An updated Momics object
      :rtype: Momics



   .. py:method:: export_sequence(output)

      Export sequence from a `.momics` repository as a `.fa` file.

      :param output: Prefix of the output fasta file
      :type output: Path

      :returns: An updated Momics object
      :rtype: Momics



   .. py:method:: export_track(track, output)

      Export a track from a `.momics` repository as a `.bw` file.

      :param track: Which track to remove
      :type track: str
      :param output: Prefix of the output bigwig file
      :type output: Path

      :returns: An updated Momics object
      :rtype: Momics



   .. py:method:: features(label = None)

      Extract table of ingested features sets.

      :returns: pd.DataFrame: A data frame listing one ingested feature set per row
                - if `label` is not None: pr.PyRanges: A PyRanges object of the specified feature set
      :rtype: - if `label` is None



   .. py:method:: ingest_chroms(chr_lengths, genome_version = '')

      Add chromosomes (and genome) information the `.momics` repository.

      :param chr_lengths: Chromosome lengths
      :type chr_lengths: dict
      :param genome_version: Genome version (default: "").                 Defaults to "".
      :type genome_version: str, optional

      :returns: An updated Momics object
      :rtype: Momics



   .. py:method:: ingest_features(features, threads = 1, max_features = 9999, tile = 32000)

      Ingest feature sets to the `.momics` repository.

      :param features: Dictionary of feature sets already imported as a PyRanges.
      :type features: dict
      :param threads: Threads to parallelize I/O. Defaults to 1.
      :type threads: int, optional
      :param max_features: Maximum number of feature sets. Defaults to 9999.
      :type max_features: int, optional
      :param tile: Tile size. Defaults to 50000.
      :type tile: int, optional
      :param compression: Compression level. Defaults to 3.
      :type compression: int, optional

      :returns: The updated Momics object
      :rtype: Momics



   .. py:method:: ingest_sequence(fasta, threads = 1, tile = 32000)

      Ingest a fasta file into a Momics repository

      :param fasta: Path to a Fasta file containing the genome reference sequence.
      :type fasta: str
      :param threads: Threads to parallelize I/O. Defaults to 1.
      :type threads: int, optional
      :param tile: Tile size for TileDB. Defaults to 50000.
      :type tile: int, optional

      :returns: The updated Momics object
      :rtype: Momics



   .. py:method:: ingest_track(coverage, track, threads = 1)

      Ingest a coverage track provided as a dictionary to a `.momics` repository.
      This method is useful when you have already computed the coverage track and
      have it in memory.

      :param coverage: Dictionary of coverage tracks. The keys are                 chromosome names and the values are numpy arrays.
      :type coverage: dict
      :param track: Label to store the track under.
      :type track: str
      :param threads: Threads to parallelize I/O. Defaults to 1.
      :type threads: int, optional

      :returns: The updated Momics object
      :rtype: Momics



   .. py:method:: ingest_tracks(bws, threads = 1, max_bws = 9999, tile = 32000)

      Ingest bigwig coverage tracks to the `.momics` repository.

      :param bws: Dictionary of bigwig files
      :type bws: dict
      :param threads: Threads to parallelize I/O. Defaults to 1.
      :type threads: int, optional
      :param max_bws: Maximum number of bigwig files. Defaults to 9999.
      :type max_bws: int, optional
      :param tile: Tile size. Defaults to 50000.
      :type tile: int, optional
      :param compression: Compression level. Defaults to 3.
      :type compression: int, optional

      :returns: The updated Momics object
      :rtype: Momics



   .. py:method:: manifest()

      Returns the manifest of the Momics repository. The manifest lists the
      configuration for all the arrays stored in the repository, including
      their schema, attributes, and metadata.



   .. py:method:: remove()

      Remove a `.momics` repository.



   .. py:method:: remove_track(track)

      Remove a track from a `.momics` repository.

      :param track: Which track to remove
      :type track: str

      :returns: An updated Momics object
      :rtype: Momics



   .. py:method:: seq(label = None)

      Extract sequence table from a `.momics` repository.

      :param label: Which chromosome to extract. Defaults to None.
      :type label: str, optional

      :returns: A data frame listing one chromosome per row,
                with first/last 10 nts.
      :rtype: pd.DataFrame



   .. py:method:: size()

      :returns: The size of the repository in bytes.
      :rtype: int



   .. py:method:: tracks(label = None)

      Extract table of ingested bigwigs.

      :returns: A data frame listing one ingested bigwig file per row
      :rtype: pd.DataFrame



   .. py:attribute:: cfg
      :value: None



   .. py:attribute:: path


.. py:data:: BULK_OVERWRITE
   :value: True


.. py:data:: TILEDB_CHUNKSIZE
   :value: 10000


.. py:data:: TILEDB_COMPRESSION
   :value: 2


.. py:data:: TILEDB_COV_FILTERS

.. py:data:: TILEDB_POSITION_FILTERS

.. py:data:: TILEDB_SEQ_FILTERS

.. py:data:: lock


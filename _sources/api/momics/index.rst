momics
======

.. py:module:: momics

.. autoapi-nested-parse::

   momics
   ~~~~~~

   Cloud-native, TileDB-based multi-omics data format.

   :author: Jacques Serizay
   :license: CC BY-NC 4.0



Subpackages
-----------

.. toctree::
   :maxdepth: 1

   /api/momics/cli/index


Submodules
----------

.. toctree::
   :maxdepth: 1

   /api/momics/chromnn/index
   /api/momics/config/index
   /api/momics/dataloader/index
   /api/momics/export/index
   /api/momics/logging/index
   /api/momics/momics/index
   /api/momics/multirangequery/index
   /api/momics/utils/index
   /api/momics/version/index


Attributes
----------

.. autoapisummary::

   momics.__format_version__
   momics.__version__


Classes
-------

.. autoapisummary::

   momics.Momics
   momics.MultiRangeQuery


Package Contents
----------------

.. py:class:: Momics(path, config = None)

   A class to manipulate `.momics` repositories.

   .. attribute:: path

      Path to a `.momics` repository.

      :type: str


   .. py:method:: bins(width, stride, cut_last_bin_out=False)

      Generate a BedTool of tiled genomic bins

      :param width: The width of each bin.
      :type width: _type_
      :param stride: The stride size for tiling.
      :type stride: _type_
      :param cut_last_bin_out: Remove the last bin of each                 chromosome. Defaults to False.
      :type cut_last_bin_out: bool, optional
      :param Remember that PyRanges are 0-based and end-exclusive.:

      :returns: pr.PyRanges: a PyRanges object of tiled genomic bins.
      :rtype: _type_



   .. py:method:: chroms()

      Extract chromosome table from a `.momics` repository.

      :returns: A data frame listing one chromosome per row
      :rtype: pd.DataFrame



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

      :param features: Dictionary of feature sets already imported with                 pyBedTools.
      :type features: dict
      :param threads: Threads to parallelize I/O. Defaults to 1.
      :type threads: int, optional
      :param max_features: Maximum number of feature sets.                 Defaults to 9999.
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



   .. py:method:: remove()

      Remove a `.momics` repository.



   .. py:method:: remove_track(track)

      Remove a track from a `.momics` repository.

      :param track: Which track to remove
      :type track: str

      :returns: An updated Momics object
      :rtype: Momics



   .. py:method:: seq()

      Extract sequence table from a `.momics` repository.

      :returns: A data frame listing one chromosome per row,
                with first/last 10 nts.
      :rtype: pd.DataFrame



   .. py:method:: tracks(label = None)

      Extract table of ingested bigwigs.

      :returns: A data frame listing one ingested bigwig file per row
      :rtype: pd.DataFrame



   .. py:attribute:: cfg


   .. py:attribute:: path


.. py:class:: MultiRangeQuery(momics, bed)

   A class to query `.momics` repositories.

   .. attribute:: momics (Momics)



      :type: a local `.momics` repository.

   .. attribute:: queries (dict)



      :type: Dict. of pr.PyRanges object.

   .. attribute:: coverage (dict)



      :type: Dictionary of coverage scores extracted from the         `.momics` repository, populated after calling `q.query_tracks()`

   .. attribute:: seq (dict)



      :type: Dictionary of sequences extracted from the `.momics`         repository, populated after calling `q.query_seq()`


   .. py:method:: query_sequence(threads = None)

      Query multiple sequence ranges from a Momics repo.

      :param threads: Number of threads for parallel query.                 Defaults to all.
      :type threads: int, optional

      :returns: An updated MultiRangeQuery object
      :rtype: MultiRangeQuery



   .. py:method:: query_tracks(threads = None, tracks = None)

      Query multiple coverage ranges from a Momics repo.

      :param threads: Number of threads for parallel query.                 Defaults to all.
      :type threads: int, optional
      :param tracks: List of tracks to query. Defaults to None,                 which queries all tracks.
      :type tracks: list, optional

      :returns: MultiRangeQuery: An updated MultiRangeQuery object
      :rtype: MultiRangeQuery



   .. py:method:: to_SeqRecord()

      Parse self.seq attribute to a SeqRecord

      :returns: `self.seq` dictionary wrangled into a SeqRecord
      :rtype: SeqRecord



   .. py:method:: to_df()

      Parse self.coverage attribute to a pd.DataFrame

      :returns: `self.coverage` dictionary wrangled into a pd.DataFrame
      :rtype: pd.DataFrame



   .. py:method:: to_json(output)

      Write the results of a multi-range query to a JSON file.

      :param output: Path to the output JSON file.
      :type output: Path



   .. py:method:: to_npz(output)

      Write the results of a multi-range query to a NPZ file.

      :param output: Path to the output NPZ file.
      :type output: Path



   .. py:attribute:: coverage
      :type:  Optional[Dict]
      :value: None



   .. py:attribute:: momics


   .. py:attribute:: ranges


   .. py:attribute:: seq
      :type:  Optional[Dict]
      :value: None



.. py:data:: __format_version__
   :value: 1


.. py:data:: __version__


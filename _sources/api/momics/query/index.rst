momics.query
============

.. py:module:: momics.query


Classes
-------

.. autoapisummary::

   momics.query.MomicsQuery


Module Contents
---------------

.. py:class:: MomicsQuery(momics, bed)

   A class to query `.momics` repositories.

   A `MomicsQuery` object can be used to query coverage and sequence data from a         `.momics` repository. It is initialized with a `Momics` object and a         `pr.PyRanges` object. The `query_tracks()` and `query_sequence()` methods         can be used to query coverage and sequence data, respectively. These methods         populate the `coverage` and `seq` attributes of the `MomicsQuery` object.

   .. attribute:: momics

      a local `.momics` repository.

      :type: Momics

   .. attribute:: queries

      `pr.PyRanges` object

      :type: pr.PyRanges

   .. attribute:: coverage

      Dictionary of coverage scores extracted from the             `.momics` repository, populated after calling `q.query_tracks()`.             Dictionary keys are the names of the tracks, and values are dictionaries             of ranges and scores. For example, `q.coverage["bw1"]["chr1:0-100"]`             returns the scores for the range `chr1:0-100` in the track `bw1`.             The scores are stored as a list of numpy arrays, where each array             corresponds to the scores for a single range.

      :type: dict

   .. attribute:: seq

      Dictionary of sequences extracted from the `.momics`             repository, populated after calling `q.query_seq()`. The structure of
      the `seq` dictionary is similar to that of the `coverage` dictionary,             but the values are strings of nucleotide sequences instead of lists of             numpy arrays. For example, `q.seq["nucleotide"]["chr1:0-100"]` returns             the nucleotide sequence for the range `chr1:0-100`. The sequences are             stored as strings, where each string corresponds to the sequence for a             single range.

      :type: dict

   :returns: A `MomicsQuery` object.
   :rtype: MomicsQuery


   .. py:method:: query_sequence(threads = None, silent = True)

      Query multiple sequence ranges from a Momics repo.

      :param threads: Number of threads for parallel query. Defaults to all.
      :type threads: int, optional
      :param silent: Whether to suppress info messages.
      :type silent: bool, optional

      :returns: An updated MomicsQuery object
      :rtype: MomicsQuery



   .. py:method:: query_tracks(threads = None, tracks = None, silent = True)

      Query multiple coverage ranges from a Momics repo.

      :param threads: Number of threads for parallel query.                 Defaults to all.
      :type threads: int, optional
      :param tracks: List of tracks to query. Defaults to None,                 which queries all tracks.
      :type tracks: list, optional
      :param silent: Whether to suppress info messages.
      :type silent: bool, optional

      :returns: MomicsQuery: An updated MomicsQuery object
      :rtype: MomicsQuery



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
      :type:  Optional[dict]
      :value: None



   .. py:attribute:: momics


   .. py:attribute:: ranges


   .. py:attribute:: seq
      :type:  Optional[dict]
      :value: None




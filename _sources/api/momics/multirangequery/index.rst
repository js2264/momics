momics.multirangequery
======================

.. py:module:: momics.multirangequery


Classes
-------

.. autoapisummary::

   momics.multirangequery.MultiRangeQuery


Module Contents
---------------

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




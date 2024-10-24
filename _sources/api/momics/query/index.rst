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

   .. attribute:: momics

      a local `.momics` repository.

      :type: Momics

   .. attribute:: queries

      `pr.PyRanges` object

      :type: pr.PyRanges

   .. attribute:: coverage

      Dictionary of coverage scores extracted from the             `.momics` repository, populated after calling `q.query_tracks()`

      :type: dict

   .. attribute:: seq

      Dictionary of sequences extracted from the `.momics`             repository, populated after calling `q.query_seq()`

      :type: dict


   .. py:method:: pileup(type = 'mean', prefix = None)

      Aggregate query coverages into genome-wide dictionary(ies).
      If the `coverage` attribute has not been populated yet, the :func:`query_tracks()` method will be called.
      The coverage over each range is aggregated across all tracks. In the case of
      overlapping ranges, the coverage is averaged.

      Each value of the output dictionary is a dictionary itself, with the keys being the chromosome names
      and the values being the coverage score, averaged for overlapping ranges.

      :param prefix: Prefix to the output `.bw` files to create.
                     If provided, queried coverage will be saved for each track in a file
                     named `<prefix>_<track_label>.bw`.
      :type prefix: str, optional

      :returns: A dictionary of genome-wide coverage scores, for each track. If
                the queried ranges overlap, the coverage is averaged.
                Note that if the output argument is provided, the results will be
                saved to a `.bw` file and this function will return the Path to the file.

      .. seealso:: :func:`MomicsQuery.query_tracks()`

      .. rubric:: Examples

      >>> mom = momics.momics.Momics('path/to/momics')
      >>> windows = pr.PyRanges(
      ...     chromosomes = ["I", "I", "I", "I"],
      ...     starts = [0, 5, 10, 20],
      ...     ends = [30, 30, 30, 30],
      ... )
      >>> q = MomicsQuery(mom, windows)
      >>> q.query_tracks()
      >>> q.pileup()



   .. py:method:: query_sequence(threads = None)

      Query multiple sequence ranges from a Momics repo.

      :param threads: Number of threads for parallel query.                 Defaults to all.
      :type threads: int, optional

      :returns: An updated MomicsQuery object
      :rtype: MomicsQuery



   .. py:method:: query_tracks(threads = None, tracks = None)

      Query multiple coverage ranges from a Momics repo.

      :param threads: Number of threads for parallel query.                 Defaults to all.
      :type threads: int, optional
      :param tracks: List of tracks to query. Defaults to None,                 which queries all tracks.
      :type tracks: list, optional

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
      :type:  Optional[Dict]
      :value: None



   .. py:attribute:: momics


   .. py:attribute:: ranges


   .. py:attribute:: seq
      :type:  Optional[Dict]
      :value: None




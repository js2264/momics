momics.aggregate
================

.. py:module:: momics.aggregate


Functions
---------

.. autoapisummary::

   momics.aggregate.aggregate


Module Contents
---------------

.. py:function:: aggregate(cov, ranges, chrom_sizes, type = 'mean', prefix = None)

   Aggregate query coverage outputs into genome-wide dictionary(ies).
   The coverage over each range is aggregated across all tracks. In the case of
   overlapping ranges, the coverage is averaged.

   Each value of the output dictionary is a dictionary itself, with the keys being the chromosome names
   and the values being the coverage score, averaged for overlapping ranges.

   :param cov: A dictionary of coverage scores, for each track. This is generally the output of
               :func:`MomicsQuery.query_tracks().coverage`.
   :type cov: dict
   :param ranges: A PyRanges object containing the ranges queried.
   :type ranges: PyRanges
   :param chrom_sizes: A dictionary of chromosome sizes.
   :type chrom_sizes: dict
   :param type: The type of aggregation to perform. Can be either "mean" or "sum".
   :param prefix: Prefix to the output `.bw` files to create.
                  If provided, queried coverage will be saved for each track in a file
                  named `<prefix>_<track_label>.bw`.
   :type prefix: str, optional

   :returns: A dictionary of genome-wide coverage scores, for each track. If
             the queried ranges overlap, the coverage is averaged/summed.
             Note that if the output argument is provided, the results for each track will be
             saved to a `<prefix>_<track_label>.bw` file.

   .. seealso:: :func:`MomicsQuery.query_tracks()`

   .. rubric:: Examples

   >>> mom = momics.momics.Momics('path/to/momics')
   >>> windows = pr.PyRanges(
   ...     chromosomes = ["I", "I", "I", "I"],
   ...     starts = [0, 5, 10, 20],
   ...     ends = [30, 30, 30, 30],
   ... )
   >>> cov = MomicsQuery(mom, windows).coverage
   >>> aggregate(cov, windows, {"I": 30})



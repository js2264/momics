momics.utils
============

.. py:module:: momics.utils


Functions
---------

.. autoapisummary::

   momics.utils.dict_to_bigwig
   momics.utils.get_chr_lengths
   momics.utils.parse_ucsc_coordinates
   momics.utils.pyranges_to_bw
   momics.utils.split_ranges


Module Contents
---------------

.. py:function:: dict_to_bigwig(bw_dict, output)

   Write a dictionary of coverages to a bigwig file.
   The dictionary should have chromosome names as keys and per-base coverage as values.

   :param bw_dict: Dictionary of chromosome coverages
   :type bw_dict: dict
   :param output: Path to output bigwig file
   :type output: Path

   :returns: Path to the output bigwig file

   .. rubric:: Examples

   >>> bw_dict = {'chr1': np.random.rand(1000), 'chr2': np.random.rand(2000)}
   >>> dict_to_bigwig(bw_dict, 'output.bw')


.. py:function:: get_chr_lengths(bw)

   A simple wrapper around pyBigWig to get chromosome lengths from a bigwig file.

   :param bw: path to a bigwig file
   :type bw: Path

   :returns: Dictionary of chromosome lengths
   :rtype: dict


.. py:function:: parse_ucsc_coordinates(coords)

   Parse UCSC-style coordinates as a pr.PyRanges object. The coordinates should be in the format "chrom:start-end".

   :param coords: A UCSC-style set of coordinates (e.g., "I:11-100").
   :type coords: str

   :returns: A pr.PyRanges object.
   :rtype: pr.PyRanges


.. py:function:: pyranges_to_bw(pyranges, scores, output)

   Write a PyRanges object and corresponding scores to a BigWig file.
   The PyRanges object must have the same length as the first dimension of the scores array.
   The PyRanges object must have ranges of the same width as the second dimension of the scores array.

   :param pyranges: A PyRanges object.
   :type pyranges: pr.PyRanges
   :param scores: A 2D NumPy array of scores.
   :type scores: np.ndarray
   :param output: Path to the output BigWig file.
   :type output: str

   :returns: None


.. py:function:: split_ranges(pyranges, ratio=0.8, shuffle=True)

   Split a PyRanges object into two PyRanges objects based on a ratio.
   The first PyRanges object will contain the first `ratio` proportion of the
   ranges, and the second PyRanges object will contain the remaining ranges.

   :param pyranges: A PyRanges object.
   :type pyranges: pr.PyRanges
   :param ratio: A float between 0 and 1.
   :type ratio: float

   :returns: A tuple of two PyRanges objects.
   :rtype: Tuple[pr.PyRanges, pr.PyRanges]



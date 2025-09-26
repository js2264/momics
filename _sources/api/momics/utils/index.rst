momics.utils
============

.. py:module:: momics.utils


Attributes
----------

.. autoapisummary::

   momics.utils.DEFAULT_OHD_MAPPING
   momics.utils.DEFAULT_OHE_MAPPING


Functions
---------

.. autoapisummary::

   momics.utils.dict_to_bigwig
   momics.utils.get_chr_lengths
   momics.utils.one_hot_decode
   momics.utils.one_hot_encode
   momics.utils.parse_ucsc_coordinates
   momics.utils.pyranges_to_bw
   momics.utils.scale_track
   momics.utils.split_ranges
   momics.utils.to_bw


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


.. py:function:: one_hot_decode(encoded_sequences, mapping = DEFAULT_OHD_MAPPING)

   Decode one-hot encoded sequences back to their original string representation.

   :param encoded_sequences: A one-hot encoded array of shape (num_sequences, seq_length, 4).
   :type encoded_sequences: np.ndarray

   :returns: A list of decoded DNA sequences.
   :rtype: List[str]

   .. rubric:: Examples

   >>> encoded = np.array([[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]])
   >>> one_hot_decode(encoded)
   ['ATGC']


.. py:function:: one_hot_encode(sequences, mapping=DEFAULT_OHE_MAPPING, handle_non_standard=False, dtype=np.int8)

   Efficiently one-hot encode DNA sequences.

   :param sequences: A single DNA sequence or list of DNA sequences
   :type sequences: Union[str, List[str]]
   :param handle_non_standard: If True, non-standard nucleotides (not A,T,G,C) will be
                               encoded as [0,0,0,0]. If False, will raise a KeyError.
   :type handle_non_standard: bool
   :param dtype: NumPy data type for the output array (default: np.int8 to save memory)

   :returns:

             A one-hot encoded array of shape (len(sequences), seq_length, 4) for multiple
                 sequences or (seq_length, 4) for a single sequence. With the default mapping, columns represent:
                 - Column 1: Adenine (A)
                 - Column 2: Thymine (T)
                 - Column 3: Guanine (G)
                 - Column 4: Cytosine (C)
   :rtype: np.ndarray

   .. rubric:: Examples

   >>> one_hot_encode("ATGC")
   array([[1, 0, 0, 0],
          [0, 1, 0, 0],
          [0, 0, 1, 0],
          [0, 0, 0, 1]], dtype=int8)
   >>> one_hot_encode(["ATGC", "TACG"])
   array([[[1, 0, 0, 0],
           [0, 1, 0, 0],
           [0, 0, 1, 0],
           [0, 0, 0, 1]],
          [[0, 1, 0, 0],
           [1, 0, 0, 0],
           [0, 0, 0, 1],
           [0, 0, 1, 0]]], dtype=int8)


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


.. py:function:: scale_track(cov, quartile = 99.99, blacklist = None, threshold = None)

   Scale coverage values in a dictionary of chromosome coverages to the range [0, 1].
   The scaling is done by dividing each value by the maximum value in the chromosome,
   after capping values at a specified percentile threshold.
   :param cov: A dictionary where keys are chromosome names and values are coverage arrays.
   :type cov: dict
   :param quartile: The *genome-wide* percentile to use for capping coverage values.
   :type quartile: float
   :param blacklist: A blacklist of regions to exclude when computing quartile.
   :type blacklist: Optional[pr.PyRanges]
   :param threshold: If provided, by-passes quartile calculation and uses this value for capping.
   :type threshold: Optional[float]

   :returns: A dictionary with scaled coverage values.
   :rtype: dict


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


.. py:function:: to_bw(dict, bw_names = None)

   Export multiple tracks into separate bigwig files. Internally, it wraps
   the `dict_to_bigwig` function on each value of the input dictionary.

   :param dict: Dictionary of chromosome coverages
   :type dict: dict
   :param bw_names: List of names (without `.bw` extension)
                    for the output bigwig files.
                    If None, the keys of the input dictionary will be used.
   :type bw_names: list, optional


.. py:data:: DEFAULT_OHD_MAPPING

.. py:data:: DEFAULT_OHE_MAPPING


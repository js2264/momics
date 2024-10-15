momics.utils
============

.. py:module:: momics.utils


Functions
---------

.. autoapisummary::

   momics.utils.get_chr_lengths
   momics.utils.parse_ucsc_coordinates


Module Contents
---------------

.. py:function:: get_chr_lengths(bw)

   Parse bigwig header to extract chromosome lengths

   :param bw: path to a bigwig file
   :type bw: Path

   :returns: Dictionary of chromosome lengths
   :rtype: dict


.. py:function:: parse_ucsc_coordinates(coords)

   Parse UCSC-style coordinates as a pr.PyRanges object.

   :param coords: A UCSC-style set of coordinates (e.g., "I:11-100").
   :type coords: str

   :returns: A pr.PyRanges object.
   :rtype: pr.PyRanges



momics.export
=============

.. py:module:: momics.export


Functions
---------

.. autoapisummary::

   momics.export.export_features
   momics.export.export_sequence
   momics.export.export_track


Module Contents
---------------

.. py:function:: export_features(momics, features, output)

   Export a features set from a `.momics` repository as a `.bed `file.

   :param features: Which features to remove
   :type features: str
   :param output: Prefix of the output BED file
   :type output: Path

   :returns: An updated Momics object
   :rtype: Momics


.. py:function:: export_sequence(momics, output)

   Export sequence from a `.momics` repository as a `.fa `file.

   :param output: Prefix of the output bigwig file
   :type output: Path

   :returns: An updated Momics object
   :rtype: Momics


.. py:function:: export_track(momics, track, output)

   Export a track from a `.momics` repository as a `.bw `file.

   :param track: Which track to remove
   :type track: str
   :param output: Prefix of the output bigwig file
   :type output: Path

   :returns: An updated Momics object
   :rtype: Momics



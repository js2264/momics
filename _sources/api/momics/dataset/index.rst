momics.dataset
==============

.. py:module:: momics.dataset


Classes
-------

.. autoapisummary::

   momics.dataset.MomicsDataset


Module Contents
---------------

.. py:class:: MomicsDataset(variant_tensor)

   Bases: :py:obj:`tensorflow.data.Dataset`


   This class is implemented to train deep learning models, where the
   input data (features) is a track or a sequence and the labeled data (target)
   is another track. The data loader will iterate over the ranges in batches
   and extract the features and target for each range. It is a subclass of
   `tf.data.DataSet` and can be used as a generator for a `tf.keras.Model`.

   For a more basic generator to stream a `momics` by batches of ranges,
   see `momics.streamer.MomicsStreamer`.



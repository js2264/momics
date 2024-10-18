momics.dataloader
=================

.. py:module:: momics.dataloader


Classes
-------

.. autoapisummary::

   momics.dataloader.RangeDataLoader


Module Contents
---------------

.. py:class:: RangeDataLoader(momics, ranges, features, target, target_size = None, batch_size = None, silent = False)

   Bases: :py:obj:`tensorflow.keras.utils.Sequence`


   This class is implemented to train deep learning models, where the
   input data (features) is a track or a sequence and the labeled data (target)
   is another track. The data loader will iterate over the ranges in batches
   and extract the features and target for each range.

   .. attribute:: momics (Momics)



      :type: a local `.momics` repository.

   .. attribute:: ranges (dict)



      :type: pr.PyRanges object.

   .. attribute:: features (str)



      :type: the name of the track to use for input data

   .. attribute:: target (str)



      :type: the name of the track to use for output data

   .. attribute:: target_size (int)



      :type: To which width should the target be centered


   .. py:method:: __getitem__(idx)


   .. py:method:: __len__()


   .. py:method:: __str__()


   .. py:attribute:: batch_size


   .. py:attribute:: current


   .. py:attribute:: df


   .. py:attribute:: features


   .. py:attribute:: momics


   .. py:attribute:: ranges


   .. py:attribute:: silent


   .. py:attribute:: start
      :value: 0



   .. py:attribute:: stop


   .. py:attribute:: target


   .. py:attribute:: target_size


   .. py:attribute:: tr


   .. py:attribute:: widths



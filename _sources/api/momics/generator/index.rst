momics.generator
================

.. py:module:: momics.generator


Classes
-------

.. autoapisummary::

   momics.generator.RangeGenerator


Module Contents
---------------

.. py:class:: RangeGenerator(momics, ranges, batch_size, data, label, label_center = None)

   This class is implemented to train deep learning models, where the
   input data is a track or a sequence and the label is another track.
   The generator will iterate over the ranges in batches and extract
   the data and label for each range.

   .. attribute:: momics (Momics)



      :type: a local `.momics` repository.

   .. attribute:: ranges (dict)



      :type: pr.PyRanges object.

   .. attribute:: data (str)



      :type: the name of the track to use as data

   .. attribute:: label (str)



      :type: the name of the track to use as label

   .. attribute:: label_center (int)



      :type: To which width should the label be centered


   .. py:method:: __iter__()


   .. py:method:: __next__()


   .. py:method:: __str__()


   .. py:method:: reset()


   .. py:attribute:: current


   .. py:attribute:: data


   .. py:attribute:: df


   .. py:attribute:: label


   .. py:attribute:: label_center


   .. py:attribute:: momics


   .. py:attribute:: ranges


   .. py:attribute:: start
      :value: 0



   .. py:attribute:: step


   .. py:attribute:: stop


   .. py:attribute:: tr


   .. py:attribute:: widths



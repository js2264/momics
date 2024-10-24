momics.streamer
===============

.. py:module:: momics.streamer


Classes
-------

.. autoapisummary::

   momics.streamer.MomicsStreamer


Module Contents
---------------

.. py:class:: MomicsStreamer(momics, ranges, batch_size = None, features = None, preprocess_func = None, silent = True)

   This class is implemented to efficiently query a `momics` repository by batches
   and extract any coverage data from it. The data streamer will iterate over ranges in batches
   and iteratively query a `momics`.

   For a tensorflow DataSet constructor, see `momics.dataset.MomicsDataset`.

   .. seealso:: :class:`momics.dataset.MomicsDataset`

   .. attribute:: momics

      a local `.momics` repository.

      :type: Momics

   .. attribute:: ranges

      pr.PyRanges object.

      :type: dict

   .. attribute:: batch_size

      the batch size

      :type: int

   .. attribute:: features

      list of track labels to query

      :type: list

   .. attribute:: silent

      whether to suppress info messages

      :type: bool


   .. py:method:: __iter__()


   .. py:method:: __len__()


   .. py:method:: __next__()

      Return the next batch or raise StopIteration.



   .. py:method:: batch(batch_size)

      Change the batch size for streaming data.

      :param batch_size: The new size for batches.
      :type batch_size: int



   .. py:method:: generator()

      Generator to yield batches of ranges and queried/preprocessed data.

      :Yields: *Tuple[pr.PyRanges, np.ndarray]* -- batch_ranges and preprocessed_data



   .. py:method:: query(batch_ranges)

      Query function to fetch data from a `momics` repo based on batch_ranges.

      :param batch_ranges: PyRanges object for a batch
      :type batch_ranges: pr.PyRanges

      :returns: Queried coverage/sequence data



   .. py:method:: reset()

      Reset the iterator to allow re-iteration.



   .. py:attribute:: batch_index
      :value: 0



   .. py:attribute:: batch_size


   .. py:attribute:: features


   .. py:attribute:: momics


   .. py:attribute:: num_batches


   .. py:attribute:: preprocess_func


   .. py:attribute:: ranges


   .. py:attribute:: silent



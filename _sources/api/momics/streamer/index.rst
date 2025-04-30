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

   This class is implemented to efficiently query a `momics` repository for a set of features by batch,
   and extract data from it. The data streamer will iterate over ranges in batches
   and iteratively query a `momics`. The streamer can also be used to preprocess the data before returning it.

   When iterating over the streamer by batch, each iteration will return a dictionary with the queried data.
   The keys of the dictionary are the queried features and the values are numpy arrays
   containing the data over the corresponding ranges.

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

   .. rubric:: Example

   >>> from momics import streamer as mms
   >>> from momics import momics as mmm
   >>> from pyranges import PyRanges
   >>>
   >>> # Create a MomicsStreamer object
   >>> repo = mmm.Momics("yeast.momics")
   >>> ranges = repo.bins(1000)
   >>> stream = mms.MomicsStreamer(repo, ranges, batch_size=1000, features=["nucleotide", "atac"])
   >>>
   >>> # Iterate over the streamer
   >>> for batch in stream:
   >>>     # Process the batch
   >>>     print(batch)


   .. py:method:: __iter__()


   .. py:method:: __len__()


   .. py:method:: __next__()

      Return the next batch or raise StopIteration.



   .. py:method:: batch(batch_size)

      Change the batch size for streaming data.

      :param batch_size: The new size for batches.
      :type batch_size: int



   .. py:method:: generator()

      Generator to yield batches of data.



   .. py:method:: query(batch_ranges)

      Query function to fetch data from a `momics` repo based on batch_ranges.

      :param batch_ranges: PyRanges object for a batch
      :type batch_ranges: pr.PyRanges

      :returns: Queried coverage/sequence data



   .. py:attribute:: batch_index
      :value: 0



   .. py:attribute:: batch_size
      :value: None



   .. py:attribute:: features
      :value: None



   .. py:attribute:: momics


   .. py:attribute:: num_batches


   .. py:attribute:: preprocess_func


   .. py:attribute:: ranges


   .. py:attribute:: silent
      :value: True




momics.chromnn
==============

.. py:module:: momics.chromnn


Attributes
----------

.. autoapisummary::

   momics.chromnn.DEFAULT_CHROMNN_INPUT_LAYER
   momics.chromnn.DEFAULT_CHROMNN_OUTPUT_LAYER
   momics.chromnn.kernel_init


Classes
-------

.. autoapisummary::

   momics.chromnn.ChromNN


Module Contents
---------------

.. py:class:: ChromNN(input=DEFAULT_CHROMNN_INPUT_LAYER, output=DEFAULT_CHROMNN_OUTPUT_LAYER)

   This class implements a convolutional neural network for the prediction of
   chromatin modality from another modality. The model implements a series of
   convolutional blocks with residual connections and dropout layers.


   .. py:attribute:: kernel_init


   .. py:attribute:: model


   .. py:attribute:: x


   .. py:attribute:: x1


   .. py:attribute:: x2


   .. py:attribute:: x3


   .. py:attribute:: x4


.. py:data:: DEFAULT_CHROMNN_INPUT_LAYER

.. py:data:: DEFAULT_CHROMNN_OUTPUT_LAYER

.. py:data:: kernel_init


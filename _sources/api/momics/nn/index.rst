momics.nn
=========

.. py:module:: momics.nn


Attributes
----------

.. autoapisummary::

   momics.nn.DEFAULT_NN_INPUT_LAYER
   momics.nn.DEFAULT_NN_OUTPUT_LAYER
   momics.nn.kernel_init


Classes
-------

.. autoapisummary::

   momics.nn.Basenji
   momics.nn.ChromNN


Functions
---------

.. autoapisummary::

   momics.nn.mae_cor


Module Contents
---------------

.. py:class:: Basenji(input=DEFAULT_NN_INPUT_LAYER, output=DEFAULT_NN_OUTPUT_LAYER)

   This class is a loose adaptation of the Basenji convolutional neural network
   for the prediction of epigenomic data from DNA sequence (Kelley et al. 2018).


   .. py:attribute:: model


.. py:class:: ChromNN(inputs, outputs, filters=None, kernel_sizes=None, pool_sizes=None, dropout_rates=None, activation='relu')

   A generic neural network that can handle multiple input and output modalities.

   This network processes each input through a separate convolutional branch,
   concatenates the branches, and then splits into multiple output heads.

   :param inputs: Dictionary of input layers
   :param outputs: Dictionary of output layers
   :param filters: List of filter counts for each conv layer (default: [64, 16, 8])
   :param kernel_sizes: List of kernel sizes for each conv layer (default: [3, 8, 80])
   :param pool_sizes: List of pooling sizes for each layer (default: [2, 2, 2])
   :param dropout_rates: List of dropout rates for each layer (default: [0.2, 0.2, 0])
   :param activation: Activation function to use


   .. py:attribute:: model


.. py:function:: mae_cor(y_true, y_pred, alpha=0.5)

   Custom loss function combining MAE and correlation.

   This loss is well-suited for nucleosome positioning tasks where both the
   absolute magnitude (MAE) and pattern of peaks (correlation) are important.

   :param y_true: Ground truth values
   :param y_pred: Predicted values
   :param alpha: Weight for MAE component (1-alpha is weight for correlation)
                 Range 0-1, where 0 is pure correlation loss and 1 is pure MAE

   :returns: MAE and correlation combined loss value


.. py:data:: DEFAULT_NN_INPUT_LAYER

.. py:data:: DEFAULT_NN_OUTPUT_LAYER

.. py:data:: kernel_init


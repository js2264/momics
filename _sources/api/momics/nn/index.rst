momics.nn
=========

.. py:module:: momics.nn


Attributes
----------

.. autoapisummary::

   momics.nn.kernel_init


Classes
-------

.. autoapisummary::

   momics.nn.Basenji
   momics.nn.ChromNN
   momics.nn.Conv1DBlock
   momics.nn.DilatedConvBlock


Functions
---------

.. autoapisummary::

   momics.nn.cor
   momics.nn.loss_cor
   momics.nn.loss_mae_cor
   momics.nn.tf_rc_augmentation


Module Contents
---------------

.. py:class:: Basenji(input_size=2048, output_size=512)

   This class is a loose adaptation of the Basenji convolutional neural network
   for the prediction of epigenomic data from DNA sequence (Kelley et al. 2018).


   .. py:attribute:: model


.. py:class:: ChromNN(inputs, outputs)

   A generic neural network that can handle multiple input and output modalities.

   This network processes each input through a separate convolutional branch,
   concatenates the branches, and then splits into multiple output heads for
   true multi-modal predictions.

   :param inputs: Dictionary of input layers; all should have shape (8192, n_channels)
   :param outputs: Dictionary of output layers; all should have shape (512, 1)


   .. py:attribute:: model


.. py:class:: Conv1DBlock(n_channels, kernel_size, activation='relu', drop_out=0.0, kernel_initializer='auto', padding='same', use_batch_norm=True, pool_size=None, **kwargs)

   Bases: :py:obj:`tensorflow.keras.layers.Layer`


   Custom layer that combines Conv1D, BatchNorm, Activation, MaxPooling, and Dropout.


   .. py:method:: call(inputs, training=None)


   .. py:method:: get_config()


   .. py:attribute:: activation
      :value: 'relu'



   .. py:attribute:: activation_layer


   .. py:attribute:: conv1d


   .. py:attribute:: drop_out
      :value: 0.0



   .. py:attribute:: kernel_initializer
      :value: 'auto'



   .. py:attribute:: kernel_size


   .. py:attribute:: n_channels


   .. py:attribute:: padding
      :value: 'same'



   .. py:attribute:: pool_size
      :value: None



   .. py:attribute:: use_batch_norm
      :value: True



.. py:class:: DilatedConvBlock(n_channels, kernel_size, dilation_rate, activation='relu', drop_out=0.0, kernel_initializer='auto', padding='same', use_batch_norm=True, pool_size=None, **kwargs)

   Bases: :py:obj:`tensorflow.keras.layers.Layer`


   Custom layer that combines dilated Conv1D, BatchNorm, Activation, MaxPooling, and Dropout.


   .. py:method:: call(inputs, training=None)


   .. py:method:: get_config()


   .. py:attribute:: activation
      :value: 'relu'



   .. py:attribute:: activation_layer


   .. py:attribute:: conv1d


   .. py:attribute:: dilation_rate


   .. py:attribute:: drop_out
      :value: 0.0



   .. py:attribute:: kernel_initializer
      :value: 'auto'



   .. py:attribute:: kernel_size


   .. py:attribute:: n_channels


   .. py:attribute:: padding
      :value: 'same'



   .. py:attribute:: pool_size
      :value: None



   .. py:attribute:: use_batch_norm
      :value: True



.. py:function:: cor(y_true, y_pred)

   Returns Pearson r correlation value (-1 to 1, higher is better) for a batch.

   :param y_true: Ground truth values
   :param y_pred: Predicted values

   :returns: Pearson r correlation value


.. py:function:: loss_cor(y_true, y_pred)

   Custom loss function for correlation.

   This loss is well-suited for when the pattern of peaks (correlation) is important.

   :param y_true: Ground truth values
   :param y_pred: Predicted values

   :returns: Correlation loss value


.. py:function:: loss_mae_cor(y_true, y_pred, alpha=0.5)

   Custom loss function combining MAE and correlation.

   This loss is well-suited for nucleosome positioning tasks where both the
   absolute magnitude (MAE) and pattern of peaks (correlation) are important.

   :param y_true: Ground truth values
   :param y_pred: Predicted values
   :param alpha: Weight for MAE component (1-alpha is weight for correlation)
                 Range 0-1, where 0 is pure correlation loss and 1 is pure MAE

   :returns: MAE and correlation combined loss value


.. py:function:: tf_rc_augmentation(inputs, outputs, swapped_cols=[1, 0, 3, 2])

   Apply reverse-complement augmentation to nucleotide sequences from inputs and/or outputs.
   One of the keys in `inputs` or `outputs` should be "nucleotide" with one-hot encoded (ATGC) sequences.

   :param inputs: Dictionary of input tensors
   :param outputs: Dictionary of output tensors
   :param swapped_cols: List defining how to swap nucleotide channels for RC (default is for ATGC, ie [1, 0, 3, 2])

   :returns: Tuple of augmented inputs and outputs


.. py:data:: kernel_init


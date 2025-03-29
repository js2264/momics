import tensorflow as tf
from tensorflow.keras import layers  # type: ignore

kernel_init = tf.keras.initializers.VarianceScaling()

DEFAULT_NN_INPUT_LAYER = layers.Input(shape=(2049, 1))
DEFAULT_NN_OUTPUT_LAYER = layers.Dense(1, activation="linear")


class ChromNN:
    """
    This class implements a convolutional neural network for the prediction of
    chromatin modality from another modality. The model consists of a series of
    convolutional blocks with residual connections and dropout layers.
    """

    def __init__(self, input=DEFAULT_NN_INPUT_LAYER, output=DEFAULT_NN_OUTPUT_LAYER) -> None:
        #
        x = layers.Conv1D(32, kernel_size=5, padding="same", activation="relu", kernel_initializer=kernel_init)(input)
        x = layers.MaxPool1D(pool_size=2, padding="same")(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        x = layers.Conv1D(32, kernel_size=5, padding="same", activation="relu", kernel_initializer=kernel_init)(x)
        x = layers.MaxPool1D(pool_size=2, padding="same")(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        x = layers.Conv1D(32, kernel_size=5, padding="same", activation="relu", kernel_initializer=kernel_init)(x)
        x = layers.MaxPool1D(pool_size=2, padding="same")(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        x = layers.Conv1D(32, kernel_size=3, padding="same", activation="relu", kernel_initializer=kernel_init)(x)
        x = layers.MaxPool1D(pool_size=2, padding="same")(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        x = layers.Conv1D(32, kernel_size=1, padding="same", activation="relu", kernel_initializer=kernel_init)(x)
        x = layers.MaxPool1D(pool_size=2, padding="same")(x)
        x = layers.BatchNormalization()(x)

        x = layers.Flatten()(x)
        x = output(x)

        self.model = tf.keras.Model(input, x)


class Basenji:  # pragma: no cover
    """
    This class is a loose adaptation of the Basenji convolutional neural network
    for the prediction of epigenomic data from DNA sequence (Kelley et al. 2018).
    """

    def __init__(self, input=DEFAULT_NN_INPUT_LAYER, output=DEFAULT_NN_OUTPUT_LAYER) -> None:

        # First PooledConvLayer
        x = layers.Conv1D(64, 12, padding="same")(input)
        x = layers.ReLU()(x)
        x = layers.MaxPooling1D(4)(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # Second PooledConvLayer
        x = layers.Conv1D(64, 5, padding="same")(x)
        x = layers.ReLU()(x)
        x = layers.MaxPooling1D(2)(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # Third PooledConvLayer
        x = layers.Conv1D(64, 5, padding="same")(x)
        x = layers.ReLU()(x)
        x = layers.MaxPooling1D(2)(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # First DilatedConvLayer
        x = layers.Conv1D(32, 5, padding="same", dilation_rate=2)(x)
        x = layers.ReLU()(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # First ResidualConcatLayer
        y = layers.Conv1D(32, 5, padding="same", dilation_rate=4)(x)
        y = layers.ReLU()(y)
        y = layers.BatchNormalization()(y)
        y = layers.Dropout(0.2)(y)
        x = layers.Concatenate()([x, y])

        # Second ResidualConcatLayer
        y = layers.Conv1D(32, 5, padding="same", dilation_rate=8)(x)
        y = layers.ReLU()(y)
        y = layers.BatchNormalization()(y)
        y = layers.Dropout(0.2)(y)
        x = layers.Concatenate()([x, y])

        # Third ResidualConcatLayer
        y = layers.Conv1D(32, 5, padding="same", dilation_rate=16)(x)
        y = layers.ReLU()(y)
        y = layers.BatchNormalization()(y)
        y = layers.Dropout(0.2)(y)
        x = layers.Concatenate()([x, y])

        x = layers.Flatten()(x)
        x = output(x)

        self.model = tf.keras.Model(input, x)


class SeqCovNN:  # pragma: no cover
    """ """

    def __init__(self, seq_input, cov_input, output=DEFAULT_NN_OUTPUT_LAYER) -> None:

        # seq_input = layers.Input(shape=(seq_length, 4), name='sequence_input')
        # cov_input = layers.Input(shape=(seq_length, 1), name='cov_input')

        # Sequence processing branch (similar to Basenji)
        x_seq = layers.Conv1D(64, 12, padding="same")(seq_input)
        x_seq = layers.ReLU()(x_seq)
        x_seq = layers.MaxPooling1D(4)(x_seq)
        x_seq = layers.BatchNormalization()(x_seq)
        x_seq = layers.Dropout(0.2)(x_seq)

        x_seq = layers.Conv1D(64, 5, padding="same")(x_seq)
        x_seq = layers.ReLU()(x_seq)
        x_seq = layers.MaxPooling1D(2)(x_seq)
        x_seq = layers.BatchNormalization()(x_seq)
        x_seq = layers.Dropout(0.2)(x_seq)

        # cov data input branch

        # cov processing branch
        x_cov = layers.Conv1D(32, 12, padding="same")(cov_input)
        x_cov = layers.ReLU()(x_cov)
        x_cov = layers.MaxPooling1D(4)(x_cov)
        x_cov = layers.BatchNormalization()(x_cov)
        x_cov = layers.Dropout(0.2)(x_cov)

        x_cov = layers.Conv1D(32, 5, padding="same")(x_cov)
        x_cov = layers.ReLU()(x_cov)
        x_cov = layers.MaxPooling1D(2)(x_cov)
        x_cov = layers.BatchNormalization()(x_cov)
        x_cov = layers.Dropout(0.2)(x_cov)

        # Merge the two branches
        # At this point, both branches should have the same spatial dimensions
        merged = layers.Concatenate()([x_seq, x_cov])

        # Continue with more shared layers
        x = layers.Conv1D(64, 5, padding="same")(merged)
        x = layers.ReLU()(x)
        x = layers.MaxPooling1D(2)(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # Dilated convolutions for capturing long-range dependencies
        x = layers.Conv1D(32, 5, padding="same", dilation_rate=2)(x)
        x = layers.ReLU()(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # Residual connections with concatenation
        y = layers.Conv1D(32, 5, padding="same", dilation_rate=4)(x)
        y = layers.ReLU()(y)
        y = layers.BatchNormalization()(y)
        y = layers.Dropout(0.2)(y)
        x = layers.Concatenate()([x, y])

        y = layers.Conv1D(32, 5, padding="same", dilation_rate=8)(x)
        y = layers.ReLU()(y)
        y = layers.BatchNormalization()(y)
        y = layers.Dropout(0.2)(y)
        x = layers.Concatenate()([x, y])

        # Flatten and output
        x = layers.Flatten()(x)
        x = output(x)

        self.model = tf.keras.Model(input, x)


class NucNN:
    """
    This class implements a convolutional neural network for nucleosome prediction.
    The model consists of three convolutional layers with respectively 64, 16 and 8 kernels
    of shape (3x4), (8x64) and (80x16). A max pooling layer of size 2 and a ReLu
    activation function is applied after each of these three convolutions. Batch normalization
    and dropout of 0.2 is applied after each convolution with stride 1.

    The model takes inputs of shape (2048, 4), the last dimension representing the four
    nucleotides, and outputs a single value corresponding to the nucleosome profile
    in the middle of the input sequence.
    """

    def __init__(self, input=DEFAULT_NN_INPUT_LAYER, output=DEFAULT_NN_OUTPUT_LAYER) -> None:

        # First convolutional layer with 64 kernels of shape (3,4)
        x = layers.Conv1D(64, kernel_size=3, padding="same")(input)
        x = layers.ReLU()(x)
        x = layers.MaxPool1D(pool_size=2)(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # Second convolutional layer with 16 kernels of shape (8,64)
        x = layers.Conv1D(16, kernel_size=8, padding="same")(x)
        x = layers.ReLU()(x)
        x = layers.MaxPool1D(pool_size=2)(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # Third convolutional layer with 8 kernels of shape (80,16)
        x = layers.Conv1D(8, kernel_size=80, padding="same")(x)
        x = layers.ReLU()(x)
        x = layers.MaxPool1D(pool_size=2)(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # Flatten and output
        x = layers.Flatten()(x)
        x = output(x)

        self.model = tf.keras.Model(input, x)

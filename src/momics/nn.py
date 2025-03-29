import tensorflow as tf
from tensorflow.keras import layers  # type: ignore

kernel_init = tf.keras.initializers.VarianceScaling()

DEFAULT_NN_INPUT_LAYER = layers.Input(shape=(2049, 1))
DEFAULT_NN_OUTPUT_LAYER = layers.Dense(1, activation="linear")


class ChromNN:
    """
    This class implements a configurable convolutional neural network for the prediction of
    chromatin modality from another modality. The model consists of a series of
    convolutional blocks with residual connections and dropout layers.

    Args:
        input: Input layer
        output: Output layer
        filters: List of filter counts for each conv layer (default: [64, 16, 8])
        kernel_sizes: List of kernel sizes for each conv layer (default: [3, 8, 80])
        pool_sizes: List of pooling sizes for each layer (default: [2, 2, 2])
        dropout_rates: List of dropout rates for each layer (default: [0.2, 0.2, 0])
        activation: Activation function to use
    """

    def __init__(
        self,
        input=DEFAULT_NN_INPUT_LAYER,
        output=DEFAULT_NN_OUTPUT_LAYER,
        filters=None,
        kernel_sizes=None,
        pool_sizes=None,
        dropout_rates=None,
        activation="relu",
    ) -> None:

        if filters is None:
            filters = [64, 16, 8]
        if kernel_sizes is None:
            kernel_sizes = [3, 8, 80]
        if pool_sizes is None:
            pool_sizes = [2, 2, 2]
        if dropout_rates is None:
            dropout_rates = [0.2, 0.2, 0]
        n_layers = len(filters)

        # If dropout_rates is a single value, expand it to a list
        if isinstance(dropout_rates, (int, float)):
            dropout_rates = [dropout_rates] * n_layers

        # Ensure all parameter lists have the same length
        assert len(dropout_rates) == n_layers, "dropout_rates must have same length as filters"
        assert len(kernel_sizes) == n_layers, "kernel_sizes must have same length as filters"
        assert len(pool_sizes) == n_layers, "pool_sizes must have same length as filters"

        x = input
        for i in range(n_layers):
            x = layers.Conv1D(
                filters[i], kernel_size=kernel_sizes[i], padding="same", activation=activation, kernel_initializer=kernel_init
            )(x)
            x = layers.MaxPool1D(pool_size=pool_sizes[i], padding="same")(x)
            x = layers.BatchNormalization()(x)

            if dropout_rates[i] > 0:
                x = layers.Dropout(dropout_rates[i])(x)

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


def mae_cor(y_true, y_pred, alpha=0.5):
    """
    Custom loss function combining MAE and correlation.

    This loss is well-suited for nucleosome positioning tasks where both the
    absolute magnitude (MAE) and pattern of peaks (correlation) are important.

    Args:
        y_true: Ground truth values
        y_pred: Predicted values
        alpha: Weight for MAE component (1-alpha is weight for correlation)
               Range 0-1, where 0 is pure correlation loss and 1 is pure MAE

    Returns:
        MAE and correlation combined loss value
    """

    def correlation_coefficient(y_true, y_pred):
        y_true = tf.reshape(y_true, [tf.shape(y_true)[0], -1])
        y_pred = tf.reshape(y_pred, [tf.shape(y_pred)[0], -1])
        x_mean = tf.reduce_mean(y_true, axis=1, keepdims=True)
        y_mean = tf.reduce_mean(y_pred, axis=1, keepdims=True)
        cov_xy = tf.reduce_mean((y_true - x_mean) * (y_pred - y_mean), axis=1)
        std_x = tf.sqrt(tf.reduce_mean(tf.square(y_true - x_mean), axis=1))
        std_y = tf.sqrt(tf.reduce_mean(tf.square(y_pred - y_mean), axis=1))
        corr = cov_xy / (std_x * std_y + tf.keras.backend.epsilon())
        return tf.reduce_mean(corr)

    cor_loss = 1.0 - correlation_coefficient(y_true, y_pred)
    mae = tf.reduce_mean(tf.abs(y_true - y_pred))

    # Combine losses with weighting
    return alpha * mae + (1.0 - alpha) * cor_loss

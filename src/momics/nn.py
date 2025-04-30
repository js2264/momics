import tensorflow as tf
from tensorflow.keras import layers  # type: ignore

kernel_init = tf.keras.initializers.VarianceScaling()

DEFAULT_NN_INPUT_LAYER = layers.Input(shape=(2049, 1))
DEFAULT_NN_OUTPUT_LAYER = layers.Dense(1, activation="linear")


class ChromNN:
    """
    A generic neural network that can handle multiple input and output modalities.

    This network processes each input through a separate convolutional branch,
    concatenates the branches, and then splits into multiple output heads.

    Args:
        inputs: Dictionary of input layers
        outputs: Dictionary of output layers
        filters: List of filter counts for each conv layer (default: [64, 16, 8])
        kernel_sizes: List of kernel sizes for each conv layer (default: [3, 8, 80])
        pool_sizes: List of pooling sizes for each layer (default: [2, 2, 2])
        dropout_rates: List of dropout rates for each layer (default: [0.2, 0.2, 0])
        activation: Activation function to use
    """

    def __init__(
        self,
        inputs,
        outputs,
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

        # Process each input through its own convolutional branch
        processed_inputs = []
        input_tensors = []

        for input_layer in inputs.values():
            input_tensors.append(input_layer)
            x = input_layer

            # Apply convolutional blocks to this input
            for i in range(n_layers):
                x = layers.Conv1D(
                    filters[i], kernel_size=kernel_sizes[i], padding="same", activation=activation, kernel_initializer=kernel_init
                )(x)
                x = layers.MaxPool1D(pool_size=pool_sizes[i], padding="same")(x)
                x = layers.BatchNormalization()(x)

                if dropout_rates[i] > 0:
                    x = layers.Dropout(dropout_rates[i])(x)

            x = layers.Flatten()(x)
            processed_inputs.append(x)

        # Concatenate all processed inputs if there are multiple
        if len(processed_inputs) > 1:
            merged = layers.Concatenate()(processed_inputs)
        else:
            merged = processed_inputs[0]

        # Create separate output heads
        output_tensors = {}

        for output_name, output_layer in outputs.items():
            output_tensors[output_name] = output_layer(merged)

        # Create model with multiple inputs and outputs
        self.model = tf.keras.Model(
            inputs={name: layer for name, layer in zip(inputs.keys(), input_tensors)}, outputs=output_tensors
        )


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
        x_mean = tf.reduce_mean(y_true, axis=1, keepdims=True)
        y_mean = tf.reduce_mean(y_pred, axis=1, keepdims=True)
        cov_xy = tf.reduce_mean((y_true - x_mean) * (y_pred - y_mean), axis=1)
        std_x = tf.sqrt(tf.reduce_mean(tf.square(y_true - x_mean), axis=1) + tf.keras.backend.epsilon())
        std_y = tf.sqrt(tf.reduce_mean(tf.square(y_pred - y_mean), axis=1) + tf.keras.backend.epsilon())
        corr = cov_xy / (std_x * std_y + tf.keras.backend.epsilon())
        return tf.reduce_mean(corr)

    cor_loss = 1.0 - correlation_coefficient(y_true, y_pred)
    mae = tf.reduce_mean(tf.abs(y_true - y_pred))

    # Combine losses with weighting
    return alpha * mae + (1.0 - alpha) * cor_loss

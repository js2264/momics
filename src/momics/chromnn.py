import tensorflow as tf
from tensorflow.keras import layers  # type: ignore

kernel_init = tf.keras.initializers.VarianceScaling()

DEFAULT_CHROMNN_INPUT_LAYER = layers.Input(shape=(2048, 1))
DEFAULT_CHROMNN_OUTPUT_LAYER = layers.Dense(64, activation="linear")


class ChromNN:
    """
    This class implements a convolutional neural network for the prediction of
    chromatin modality from another modality. The model consists of a series of
    convolutional blocks with residual connections and dropout layers.
    """

    def __init__(self, input=DEFAULT_CHROMNN_INPUT_LAYER, output=DEFAULT_CHROMNN_OUTPUT_LAYER) -> None:
        kernel_init = tf.keras.initializers.VarianceScaling()

        # First convolutional block
        x = layers.Conv1D(32, kernel_size=12, padding="same", activation="relu", kernel_initializer=kernel_init)(input)
        x = layers.MaxPool1D(pool_size=8, padding="same")(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # Second convolutional block
        x = layers.Conv1D(32, kernel_size=5, padding="same", activation="relu", kernel_initializer=kernel_init)(x)
        x = layers.MaxPool1D(pool_size=4, padding="same")(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # Third convolutional block
        x = layers.Conv1D(32, kernel_size=5, padding="same", activation="relu", kernel_initializer=kernel_init)(x)
        x = layers.MaxPool1D(pool_size=4, padding="same")(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # Fourth convolutional block with residual connection
        x = layers.Conv1D(16, kernel_size=5, padding="same", activation="relu", kernel_initializer=kernel_init, dilation_rate=2)(
            x
        )
        x = layers.BatchNormalization()(x)
        x1 = layers.Dropout(0.2)(x)

        # Fifth convolutional block with residual connection
        x = x1
        x = layers.Conv1D(16, kernel_size=5, padding="same", activation="relu", kernel_initializer=kernel_init, dilation_rate=4)(
            x
        )
        x = layers.BatchNormalization()(x)
        x2 = layers.Dropout(0.2)(x)

        # Sixth convolutional block with residual connection
        x = layers.concatenate([x1, x2], axis=2)
        x = layers.Conv1D(16, kernel_size=5, padding="same", activation="relu", kernel_initializer=kernel_init, dilation_rate=8)(
            x
        )
        x = layers.BatchNormalization()(x)
        x3 = layers.Dropout(0.2)(x)

        # Seventh convolutional block with residual connection
        x = layers.concatenate([x1, x2, x3], axis=2)
        x = layers.Conv1D(16, kernel_size=5, padding="same", activation="relu", kernel_initializer=kernel_init, dilation_rate=16)(
            x
        )
        x = layers.BatchNormalization()(x)
        x4 = layers.Dropout(0.2)(x)

        # Final layers
        x = layers.concatenate([x1, x2, x3, x4], axis=2)
        x = layers.Flatten()(x)
        x = output(x)

        self.model = tf.keras.Model(input, x)

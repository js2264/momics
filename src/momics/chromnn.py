import tensorflow as tf
from tensorflow.keras import layers  # type: ignore

kernel_init = tf.keras.initializers.VarianceScaling()

DEFAULT_CHROMNN_INPUT_LAYER = layers.Input(shape=(2049, 1))
DEFAULT_CHROMNN_OUTPUT_LAYER = layers.Dense(1, activation="linear")


class ChromNN:
    """
    This class implements a convolutional neural network for the prediction of
    chromatin modality from another modality. The model consists of a series of
    convolutional blocks with residual connections and dropout layers.
    """

    def __init__(self, input=DEFAULT_CHROMNN_INPUT_LAYER, output=DEFAULT_CHROMNN_OUTPUT_LAYER) -> None:

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

        x = layers.Flatten()(x)
        x = output(x)
        x = tf.keras.layers.Reshape((1, 1))(x)

        self.model = tf.keras.Model(input, x)


class ChromNNlong:
    """
    This class implements a convolutional neural network for the prediction of
    chromatin modality from another modality. The model consists of a series of
    convolutional blocks with residual connections and dropout layers.
    """

    def __init__(self, input=DEFAULT_CHROMNN_INPUT_LAYER, output=DEFAULT_CHROMNN_OUTPUT_LAYER) -> None:

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

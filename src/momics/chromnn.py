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


class Basenji:
    """
    This class implements a convolutional neural network for the prediction of
    chromatin modality from a genomic sequence. The model consists of a series of
    convolutional blocks with residual connections and dropout layers.
    """

    def __init__(self, input=DEFAULT_CHROMNN_INPUT_LAYER, output=DEFAULT_CHROMNN_OUTPUT_LAYER) -> None:

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

        # Final Conv1D and ReLU
        x = layers.Conv1D(1, 1)(x)
        x = layers.ReLU()(x)

        # Additional Conv1D layer to ensure output shape (None, X, 1)
        x = layers.Flatten()(x)
        x = output(x)

        self.model = tf.keras.Model(input, x)


class Westbrook_basenji:
    """Derived from the BassenjiMultiNetwork2 architecture in Axel Westbrook's repository.

    It was originally implemented and trained in PyTorch, and its summary obtained as follows:

    ```
    n_tracks = 1
    architecture = "BassenjiMultiNetwork2"
    device = "cuda"
    model_file = "SCerevisiae_chromatin_NN_prediction/Trainedmodels/model_myco_nuc_pt8/model_state.pt"
    model = models.ARCHITECTURES[architecture](n_tracks).to(device)
    model.load_state_dict(torch.load(model_file))
    ```

    It was then refactored in TensorFlow using copilot.
    """

    def __init__(self, input=DEFAULT_CHROMNN_INPUT_LAYER, output=DEFAULT_CHROMNN_OUTPUT_LAYER) -> None:

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

        # Final Conv1D and ReLU
        x = layers.Conv1D(1, 1)(x)
        x = layers.ReLU()(x)

        x = output(x)

        self.model = tf.keras.Model(input, x)

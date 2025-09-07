import tensorflow as tf
from tensorflow.keras import layers  # type: ignore

kernel_init = tf.keras.initializers.VarianceScaling()


class ChromNN:
    """
    A generic neural network that can handle multiple input and output modalities.

    This network processes each input through a separate convolutional branch,
    concatenates the branches, and then splits into multiple output heads for
    true multi-modal predictions.

    Args:
        inputs: Dictionary of input layers; all should have shape (8192, n_channels)
        outputs: Dictionary of output layers; all should have shape (512, 1)
    """

    def __init__(
        self,
        inputs,
        outputs,
    ) -> None:

        input_branches: list = []

        for in_name, in_layer in inputs.items():

            if in_name == "nucleotide":
                ## Conv tower for nucleotide input
                x = Conv1DBlock(64, 3, "relu", drop_out=0.2, name=f"{in_name}_branch_conv1d_64_3")(in_layer)
                x = Conv1DBlock(64, 8, "relu", drop_out=0.2, name=f"{in_name}_branch_conv1d_64_8")(x)
                x = Conv1DBlock(64, 31, "relu", drop_out=0.2, name=f"{in_name}_branch_conv1d_64_31")(x)

                ## Dilated tower
                dilated_features = [x]
                for dilation in [2, 4, 8, 16, 32, 64, 128]:
                    x = DilatedConvBlock(32, 31, dilation, "relu", drop_out=0.2, name=f"{in_name}_branch_dconv_{dilation}")(x)
                    dilated_features.append(x)

                x = layers.Concatenate()(dilated_features)
                input_branches.append(x)

            if in_name != "nucleotide":

                ## Conv tower
                x = Conv1DBlock(32, 3, "relu", drop_out=0.1, name=f"{in_name}_branch_conv1d_32_3")(in_layer)
                x = Conv1DBlock(32, 5, "relu", drop_out=0.1, name=f"{in_name}_branch_conv1d_32_5")(x)
                x = Conv1DBlock(48, 7, "relu", drop_out=0.15, name=f"{in_name}_branch_conv1d_48_7")(x)
                x = Conv1DBlock(48, 15, "relu", drop_out=0.2, name=f"{in_name}_branch_conv1d_48_15")(x)

                # Dilated convolutions to capture multi-scale accessibility patterns
                dilated_features = [x]
                for dilation in [2, 4, 8, 16, 32, 64]:
                    x = DilatedConvBlock(24, 15, dilation, "relu", drop_out=0.2, name=f"{in_name}_branch_dconv_{dilation}")(x)
                    dilated_features.append(x)

                x = layers.Concatenate()(dilated_features)
                x = Conv1DBlock(64, 1, "relu", drop_out=0.2)(x)
                input_branches.append(x)

        # Merge input branches
        if len(input_branches) > 1:
            x = layers.Concatenate(axis=-1)(input_branches)
        else:
            x = input_branches[0]

        # Trunk
        x = layers.Conv1D(256, kernel_size=1, name="trunk_conv1d_256")(x)
        x = layers.ReLU(name="trunk_relu")(x)
        crop_size = (8192 - 512) // 2
        x = layers.Cropping1D(cropping=(crop_size, crop_size), name="trunk_crop1D")(x)

        ## Head
        output_heads = {}
        for out_name, out_layer in outputs.items():
            head_out = layers.Conv1D(128, kernel_size=1, name=f"{out_name}_head_conv1d_128")(x)
            head_out = layers.ReLU(name=f"{out_name}_head_relu")(head_out)
            head_out = layers.Conv1D(1, kernel_size=1, name=f"{out_name}_head_conv1d_1")(head_out)
            head_out = layers.Activation("softplus", name=f"{out_name}_head_softplus")(head_out)
            output_heads[out_name] = out_layer(head_out)

        # Create model with multiple inputs and outputs
        self.model = tf.keras.Model(inputs=inputs, outputs=output_heads)


class Conv1DBlock(layers.Layer):
    """
    Custom layer that combines Conv1D, BatchNorm, Activation, and Dropout.
    """

    def __init__(
        self,
        n_channels,
        kernel_size,
        activation="relu",
        drop_out=0.0,
        kernel_initializer="auto",
        padding="same",
        use_batch_norm=True,
        **kwargs,
    ):
        super(Conv1DBlock, self).__init__(**kwargs)

        self.n_channels = n_channels
        self.kernel_size = kernel_size
        self.activation = activation
        self.drop_out = drop_out
        self.kernel_initializer = kernel_initializer
        self.padding = padding
        self.use_batch_norm = use_batch_norm

        # Automatically pick a suitable initializer if not provided
        if kernel_initializer == "auto":
            if activation in ["relu", "leaky_relu", "elu", "prelu"]:
                kernel_initializer = tf.keras.initializers.HeNormal()
            else:
                kernel_initializer = tf.keras.initializers.GlorotNormal()

        # Create layers
        self.conv1d = layers.Conv1D(n_channels, kernel_size, padding=padding, kernel_initializer=kernel_initializer)
        if use_batch_norm:
            self.batch_norm = layers.BatchNormalization()

        self.activation_layer = layers.Activation(activation)
        if drop_out > 0:
            self.dropout = layers.Dropout(drop_out)

    def call(self, inputs, training=None):
        x = self.conv1d(inputs)
        if self.use_batch_norm:
            x = self.batch_norm(x, training=training)

        x = self.activation_layer(x)
        if self.drop_out > 0:
            x = self.dropout(x, training=training)

        return x

    def get_config(self):
        config = super().get_config()
        config.update(
            {
                "n_channels": self.n_channels,
                "kernel_size": self.kernel_size,
                "activation": self.activation,
                "drop_out": self.drop_out,
                "kernel_initializer": self.kernel_initializer,
                "padding": self.padding,
                "use_batch_norm": self.use_batch_norm,
            }
        )
        return config


class DilatedConvBlock(layers.Layer):
    """
    Custom layer that combines dilated Conv1D, BatchNorm, Activation, and Dropout.
    """

    def __init__(
        self,
        n_channels,
        kernel_size,
        dilation_rate,
        activation="relu",
        drop_out=0.0,
        kernel_initializer="auto",
        padding="same",
        use_batch_norm=True,
        **kwargs,
    ):
        super(DilatedConvBlock, self).__init__(**kwargs)

        self.n_channels = n_channels
        self.kernel_size = kernel_size
        self.dilation_rate = dilation_rate
        self.activation = activation
        self.drop_out = drop_out
        self.kernel_initializer = kernel_initializer
        self.padding = padding
        self.use_batch_norm = use_batch_norm

        # Automatically pick a suitable initializer if not provided
        if kernel_initializer == "auto":
            if activation in ["relu", "leaky_relu", "elu", "prelu"]:
                kernel_initializer = tf.keras.initializers.HeNormal()
            else:
                kernel_initializer = tf.keras.initializers.GlorotNormal()

        # Create layers
        self.conv1d = layers.Conv1D(
            n_channels, kernel_size, padding=padding, dilation_rate=dilation_rate, kernel_initializer=kernel_initializer
        )
        if use_batch_norm:
            self.batch_norm = layers.BatchNormalization()

        self.activation_layer = layers.Activation(activation)
        if drop_out > 0:
            self.dropout = layers.Dropout(drop_out)

    def call(self, inputs, training=None):
        x = self.conv1d(inputs)
        if self.use_batch_norm:
            x = self.batch_norm(x, training=training)

        x = self.activation_layer(x)
        if self.drop_out > 0:
            x = self.dropout(x, training=training)

        return x

    def get_config(self):
        config = super().get_config()
        config.update(
            {
                "n_channels": self.n_channels,
                "kernel_size": self.kernel_size,
                "dilation_rate": self.dilation_rate,
                "activation": self.activation,
                "drop_out": self.drop_out,
                "kernel_initializer": self.kernel_initializer,
                "padding": self.padding,
                "use_batch_norm": self.use_batch_norm,
            }
        )
        return config


def loss_mae_cor(y_true, y_pred, alpha=0.5):
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
    cor_loss = 1.0 - cor(y_true, y_pred)
    mae = tf.reduce_mean(tf.abs(y_true - y_pred))
    return alpha * mae + (1.0 - alpha) * cor_loss


def loss_cor(y_true, y_pred):
    """
    Custom loss function for correlation.

    This loss is well-suited for when the pattern of peaks (correlation) is important.

    Args:
        y_true: Ground truth values
        y_pred: Predicted values

    Returns:
        Correlation loss value
    """
    cor_loss = 1.0 - cor(y_true, y_pred)
    return cor_loss


def cor(y_true, y_pred):
    """
    Returns Pearson r correlation value (-1 to 1, higher is better) for a batch.

    Args:
        y_true: Ground truth values
        y_pred: Predicted values
    Returns:
        Pearson r correlation value
    """
    if y_true.shape[-1] == 1:
        y_true = tf.squeeze(y_true, axis=-1)
    if y_pred.shape[-1] == 1:
        y_pred = tf.squeeze(y_pred, axis=-1)
    noise = tf.random.normal(shape=tf.shape(y_pred), mean=0.0, stddev=1e-6)
    y_pred = y_pred + noise
    x_mean = tf.reduce_mean(y_true, axis=1, keepdims=True)
    y_mean = tf.reduce_mean(y_pred, axis=1, keepdims=True)
    cov_xy = tf.reduce_mean((y_true - x_mean) * (y_pred - y_mean), axis=1)
    std_x = tf.sqrt(tf.reduce_mean(tf.square(y_true - x_mean), axis=1) + tf.keras.backend.epsilon())
    std_y = tf.sqrt(tf.reduce_mean(tf.square(y_pred - y_mean), axis=1) + tf.keras.backend.epsilon())
    corr = cov_xy / (std_x * std_y + tf.keras.backend.epsilon())
    return tf.reduce_mean(corr)


def tf_rc_augmentation(inputs, outputs, swapped_cols=[1, 0, 3, 2]):
    """
    Apply reverse-complement augmentation to nucleotide sequences from inputs and/or outputs.
    One of the keys in `inputs` or `outputs` should be "nucleotide" with one-hot encoded (ATGC) sequences.

    Args:
        inputs: Dictionary of input tensors
        outputs: Dictionary of output tensors
        swapped_cols: List defining how to swap nucleotide channels for RC (default is for ATGC, ie [1, 0, 3, 2])

    Returns:
        Tuple of augmented inputs and outputs
    """
    inp = next(iter(inputs.keys()))
    batch_size = tf.shape(inputs[inp])[0]
    apply_rc = tf.random.uniform([batch_size]) < 0.5
    apply_rc = tf.reshape(apply_rc, [-1, 1, 1])

    # INPUT: RC sequences and Reverse tracks
    inputs_augmented = {}
    for in_name, in_data in inputs.items():
        if in_name == "nucleotide":
            nucleotide_r = tf.reverse(in_data, axis=[1])
            nucleotide_rc = tf.gather(nucleotide_r, swapped_cols, axis=-1)
            nucleotide_augmented = tf.where(apply_rc, nucleotide_rc, in_data)
            inputs_augmented[in_name] = nucleotide_augmented
        else:
            track_r = tf.reverse(in_data, axis=[1])
            track_augmented = tf.where(apply_rc, track_r, in_data)
            inputs_augmented[in_name] = track_augmented

    # OUTPUT: Reverse tracks
    outputs_augmented = {}
    for out_name, out_data in outputs.items():
        if out_name == "nucleotide":
            nucleotide_r = tf.reverse(out_data, axis=[1])
            nucleotide_rc = tf.gather(nucleotide_r, swapped_cols, axis=-1)
            nucleotide_augmented = tf.where(apply_rc, nucleotide_rc, out_data)
            outputs_augmented[out_name] = nucleotide_augmented
        else:
            track_r = tf.reverse(out_data, axis=[1])
            track_augmented = tf.where(apply_rc, track_r, out_data)
            outputs_augmented[out_name] = track_augmented

    return inputs_augmented, outputs_augmented

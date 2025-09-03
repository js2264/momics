import tensorflow as tf
from tensorflow.keras import layers  # type: ignore

kernel_init = tf.keras.initializers.VarianceScaling()
k_init = tf.keras.initializers.VarianceScaling()

DEFAULT_NN_INPUT_LAYER = layers.Input(shape=(2048, 1))
DEFAULT_NN_OUTPUT_LAYER = layers.Dense(1, activation="linear")


## ---------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------- ##
##                                   NETWORKS                                   ##
## ---------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------- ##


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

    def __init__(self, input_size=2048, output_size=512) -> None:

        # First PooledConvLayer
        input = layers.Input(shape=(input_size, 4))
        x = layers.Conv1D(64, 15, padding="same")(input)
        x = layers.ReLU()(x)
        x = layers.MaxPooling1D(4)(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # Second PooledConvLayer
        x = layers.Conv1D(64, 5, padding="same")(x)
        x = layers.ReLU()(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # Third PooledConvLayer
        x = layers.Conv1D(64, 5, padding="same")(x)
        x = layers.ReLU()(x)
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

        # Final layers
        x = layers.Conv1D(1, 1, padding="same")(x)
        P = x.shape[1] // output_size
        if P > 1:
            x = layers.AveragePooling1D(pool_size=P)(x)
        x = layers.Reshape((output_size,))(x)
        output = layers.Dense(output_size, activation="linear")(x)

        self.model = tf.keras.Model(input, output)


class BasenjiMulti:  # pragma: no cover
    def __init__(self, inputs, outputs, features_size, target_size) -> None:

        # init parameters
        current_length = features_size
        x = inputs["nucleotide"]

        # First PooledConvLayer
        x = layers.Conv1D(64, 15, padding="same")(x)
        x = layers.ReLU()(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # Second PooledConvLayer
        x = layers.Conv1D(64, 5, padding="same")(x)
        x = layers.ReLU()(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # Third PooledConvLayer
        x = layers.Conv1D(64, 5, padding="same")(x)
        x = layers.ReLU()(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # DilatedConvLayer
        x = layers.Conv1D(32, 5, padding="same", dilation_rate=2)(x)
        x = layers.ReLU()(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        # First DilatedConvLayer with ResidualConcat
        y = layers.Conv1D(32, 5, padding="same", dilation_rate=4)(x)
        y = layers.ReLU()(y)
        y = layers.BatchNormalization()(y)
        y = layers.Dropout(0.2)(y)
        x = layers.Concatenate()([x, y])

        # Second DilatedConvLayer with ResidualConcat
        y = layers.Conv1D(32, 5, padding="same", dilation_rate=8)(x)
        y = layers.ReLU()(y)
        y = layers.BatchNormalization()(y)
        y = layers.Dropout(0.2)(y)
        x = layers.Concatenate()([x, y])

        # Third DilatedConvLayer with ResidualConcat
        y = layers.Conv1D(32, 5, padding="same", dilation_rate=16)(x)
        y = layers.ReLU()(y)
        y = layers.BatchNormalization()(y)
        y = layers.Dropout(0.2)(y)
        x = layers.Concatenate()([x, y])

        # Trunk
        x = layers.Conv1D(256, kernel_size=1)(x)
        x = layers.ReLU()(x)

        # Crop
        if current_length > target_size:
            crop_size = (current_length - target_size) // 2
            x = layers.Cropping1D(cropping=(crop_size, crop_size))(x)
        elif current_length < target_size:
            x = layers.UpSampling1D(size=target_size // current_length)(x)

        # Separate head for each track
        output_heads = {}
        for out_name, out_layer in outputs.items():
            head_out = layers.Conv1D(128, kernel_size=1)(x)
            head_out = layers.ReLU()(head_out)
            head_out = layers.Conv1D(1, kernel_size=1)(head_out)
            head_out = layers.Activation("softplus")(head_out)
            output_heads[out_name] = out_layer(head_out)
            track_out = layers.Conv1D(1, 1, padding="same")(x)
            track_out = layers.Reshape((target_size,))(track_out)
            output_heads[out_name] = out_layer(track_out)

        self.model = tf.keras.Model(inputs=inputs, outputs=output_heads)


class Enformer:  # pragma: no cover
    def __init__(
        self,
        inputs,
        outputs,
        features_size,
        target_size,
        conv_channels=[64, 128, 256],
        conv_kernels=[15, 5, 5],
        conv_dropout_rates=[0.0, 0.1, 0.2],
        conv_pools=[2, 2, 2],
        dilations=[2, 4, 8],
        dilated_kernel=[3, 3, 3],
        transformer_blocks=8,
        transformer_heads=8,
        transformer_dropout=0.1,
        attention_dropout=0.1,
        final_conv_filters=256,
    ):

        # Put these in init parameters
        current_length = features_size
        x = inputs["nucleotide"]

        # Conv tower
        for i, _ in enumerate(conv_channels):
            x = layers.Conv1D(conv_channels[i], kernel_size=conv_kernels[i], padding="same")(x)
            x = layers.BatchNormalization()(x)
            x = layers.ReLU()(x)
            if conv_dropout_rates[i] > 0:
                x = layers.Dropout(conv_dropout_rates[i])(x)
            if conv_pools[i] > 0:
                current_length = current_length // conv_pools[i]
                x = layers.MaxPooling1D(pool_size=conv_pools[i])(x)

        # Dilated convolutions with residual connections
        for j, dilation in enumerate(dilations):
            conv = layers.Conv1D(conv_channels[-1], kernel_size=dilated_kernel[j], dilation_rate=dilation, padding="same")(x)
            conv = layers.BatchNormalization()(conv)
            conv = layers.ReLU()(conv)
            x = layers.Add()([x, conv])

        # Transformer blocks
        # x = tf.keras.layers.Lambda(lambda t: t, output_shape=convtower_out_shape[1:])(x)
        for _ in range(transformer_blocks):
            transformer = TransformerBlock(
                num_heads=transformer_heads,
                embed_dim=conv_channels[-1],
                ff_dim=conv_channels[-1] * 4,
                dropout=transformer_dropout,
                attention_dropout=attention_dropout,
            )
            x = transformer(x)

        # Final conv reduction
        x = layers.Conv1D(final_conv_filters, kernel_size=1)(x)
        x = layers.ReLU()(x)

        # Cropping to desired output length
        if current_length > target_size:
            crop_size = (current_length - target_size) // 2
            x = layers.Cropping1D(cropping=(crop_size, crop_size))(x)
        elif current_length < target_size:
            x = layers.UpSampling1D(size=target_size // current_length)(x)

        # Multi-head outputs
        output_heads = {}
        for out_name, out_layer in outputs.items():
            head_out = layers.Conv1D(128, kernel_size=1)(x)
            head_out = layers.ReLU()(head_out)
            head_out = layers.Conv1D(1, kernel_size=1)(head_out)
            head_out = layers.Activation("softplus")(head_out)
            output_heads[out_name] = out_layer(head_out)

        self.model = tf.keras.Model(inputs=inputs, outputs=output_heads)


class OriginalEnformer:  # pragma: no cover
    def __init__(
        self,
        inputs,
        outputs,
        features_size,
        target_size,
        conv_channels=[64, 128, 256],
        conv_pools=[2, 2, 2],
        dilations=[2, 4, 8],
        dilated_kernel=[3, 3, 3],
        transformer_blocks=8,
        transformer_heads=8,
        transformer_dropout=0.1,
        attention_dropout=0.1,
        final_conv_filters=256,
    ):

        # Put these in init parameters
        current_length = features_size
        x = inputs["nucleotide"]

        # Stem
        x = ConvBlock(features_size // 2, 15, 1, dilation=0)(x)
        x = RConvBlock(features_size // 2, 1, dilation=0)(x)
        x = layers.MaxPooling1D(pool_size=2)(x)

        # Conv tower
        for i, _ in enumerate(conv_channels):
            x = ConvBlock(conv_channels[i], 5, 1, dilation=dilations[i])(x)
            x = RConvBlock(conv_channels[i], 1, dilation=dilations[i])(x)
            # if conv_dropout_rates[i] > 0:
            #     x = layers.Dropout(conv_dropout_rates[i])(x)
            if conv_pools[i] > 0:
                current_length = current_length // conv_pools[i]
                x = layers.MaxPooling1D(pool_size=conv_pools[i])(x)

        # Dilated convolutions with residual connections
        for j, dilation in enumerate(dilations):
            conv = layers.Conv1D(conv_channels[-1], kernel_size=dilated_kernel[j], dilation_rate=dilation, padding="same")(x)
            conv = layers.BatchNormalization()(conv)
            conv = layers.ReLU()(conv)
            x = layers.Add()([x, conv])

        # Transformer blocks
        # x = tf.keras.layers.Lambda(lambda t: t, output_shape=convtower_out_shape[1:])(x)
        for _ in range(transformer_blocks):
            transformer = TransformerBlock(
                num_heads=transformer_heads,
                embed_dim=conv_channels[-1],
                ff_dim=conv_channels[-1] * 4,
                dropout=transformer_dropout,
                attention_dropout=attention_dropout,
            )
            x = transformer(x)

        # Final conv reduction
        x = layers.Conv1D(final_conv_filters, kernel_size=1)(x)
        x = layers.ReLU()(x)

        # Cropping to desired output length
        if current_length > target_size:
            crop_size = (current_length - target_size) // 2
            x = layers.Cropping1D(cropping=(crop_size, crop_size))(x)
        elif current_length < target_size:
            x = layers.UpSampling1D(size=target_size // current_length)(x)

        # Multi-head outputs
        output_heads = {}
        for out_name, out_layer in outputs.items():
            head_out = layers.Conv1D(128, kernel_size=1)(x)
            head_out = layers.ReLU()(head_out)
            head_out = layers.Conv1D(1, kernel_size=1)(head_out)
            head_out = layers.Activation("softplus")(head_out)
            output_heads[out_name] = out_layer(head_out)

        self.model = tf.keras.Model(inputs=inputs, outputs=output_heads)


class Scerevisiae_Scc1:
    def __init__(
        self,
        input={"nucleotide": layers.Input(shape=(32768, 4))},
        output={"SCC1": layers.Dense(128, activation="relu", name="SCC1")},
    ):

        k_init = tf.keras.initializers.VarianceScaling()
        x = input["nucleotide"]

        x = layers.Conv1D(32, kernel_size=12, padding="same", activation="relu", kernel_initializer=k_init)(x)
        x = layers.MaxPool1D(pool_size=8, padding="same")(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        x = layers.Conv1D(32, kernel_size=5, padding="same", activation="relu", kernel_initializer=k_init)(x)
        x = layers.MaxPool1D(pool_size=4, padding="same")(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        x = layers.Conv1D(32, kernel_size=5, padding="same", activation="relu", kernel_initializer=k_init)(x)
        x = layers.MaxPool1D(pool_size=4, padding="same")(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)

        x = layers.Conv1D(16, kernel_size=5, padding="same", activation="relu", kernel_initializer=k_init, dilation_rate=2)(x)
        x = layers.BatchNormalization()(x)
        x1 = layers.Dropout(0.2)(x)

        x = x1
        x = layers.Conv1D(16, kernel_size=5, padding="same", activation="relu", kernel_initializer=k_init, dilation_rate=4)(x)
        x = layers.BatchNormalization()(x)
        x2 = layers.Dropout(0.2)(x)

        x = layers.concatenate([x1, x2], axis=2)
        x = layers.Conv1D(16, kernel_size=5, padding="same", activation="relu", kernel_initializer=k_init, dilation_rate=8)(x)
        x = layers.BatchNormalization()(x)
        x3 = layers.Dropout(0.2)(x)

        x = layers.concatenate([x1, x2, x3], axis=2)
        x = layers.Conv1D(16, kernel_size=5, padding="same", activation="relu", kernel_initializer=k_init, dilation_rate=16)(x)
        x = layers.BatchNormalization()(x)
        x4 = layers.Dropout(0.2)(x)

        x = layers.concatenate([x1, x2, x3, x4], axis=2)
        output = layers.Conv1D(1, kernel_size=1, padding="same", activation="relu", kernel_initializer=k_init)(x)

        model = tf.keras.Model(input, output)
        self.model = model


class Scerevisiae_MNase:
    def __init__(self, features_size):

        k_init = tf.keras.initializers.VarianceScaling()
        input = layers.Input(shape=(features_size, 4))
        x = input

        x = layers.Conv1D(64, kernel_size=3, padding="same", activation="relu", kernel_initializer=k_init)(x)
        x = layers.MaxPool1D(2, padding="same")(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)
        x = layers.Conv1D(16, kernel_size=8, padding="same", activation="relu", kernel_initializer=k_init)(x)
        x = layers.MaxPool1D(2, padding="same")(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.2)(x)
        x = layers.Conv1D(8, kernel_size=80, padding="same", activation="relu", kernel_initializer=k_init)(x)
        x = layers.MaxPool1D(2, padding="same")(x)
        x = layers.BatchNormalization()(x)
        x = layers.Flatten()(x)
        x = layers.Dense(1, activation="relu")(x)
        output = layers.Reshape((1, 1))(x)

        model = tf.keras.Model(input, output)
        self.model = model


## ---------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------- ##
##                                   COMPONENTS                                 ##
## ---------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------- ##


class MHABlock(layers.Layer):
    def __init__(self, num_heads, embed_dim, attention_dropout=0.1):
        super().__init__()
        self.layernorm1 = layers.LayerNormalization(epsilon=1e-6)
        self.att = layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim // num_heads, dropout=attention_dropout)
        self.drop = layers.Dropout(attention_dropout)

    def call(self, inputs):
        norm = self.layernorm1(inputs)
        attn_output = self.att(norm, norm)
        drop = self.drop(attn_output)
        return inputs + drop


class TransformerBlock(layers.Layer):
    def __init__(self, num_heads, embed_dim, ff_dim, dropout=0.1, attention_dropout=0.1):
        super().__init__()
        self.att = layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim // num_heads, dropout=attention_dropout)
        self.ffn = tf.keras.Sequential(
            [layers.Dense(ff_dim, activation="relu"), layers.Dropout(dropout), layers.Dense(embed_dim), layers.Dropout(dropout)]
        )
        self.layernorm1 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = layers.LayerNormalization(epsilon=1e-6)

    def call(self, inputs):
        attn_output = self.att(inputs, inputs)
        out1 = self.layernorm1(inputs + attn_output)
        ffn_output = self.ffn(out1)
        return self.layernorm2(out1 + ffn_output)


class ConvBlock(layers.Layer):
    def __init__(self, filters, kernel_size, dilation, padding="same"):
        super().__init__()
        self.bn = layers.BatchNormalization()
        self.act = layers.GeLU()
        self.conv = layers.Conv1D(filters, kernel_size, padding=padding)

    def call(self, inputs):
        x = self.bn(inputs)
        x = self.act(x)
        x = self.conv(x)
        return x


class RConvBlock(layers.Layer):
    def __init__(self, filters, kernel_size, dilation, padding="same"):
        super().__init__()
        self.bn = layers.BatchNormalization()
        self.act = layers.GeLU()
        self.conv = layers.Conv1D(filters, kernel_size, padding=padding)

    def call(self, inputs):
        x = self.bn(inputs)
        x = self.act(x)
        x = self.conv(x)
        x = layers.Add()([x, inputs])
        return x


## ---------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------- ##
##                                   UTILS                                      ##
## ---------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------- ##


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

import pytest

import numpy as np
import tensorflow as tf
from momics import momics
from momics import query
from momics import nn
from tensorflow.keras import layers  # type: ignore
from momics.aggregate import aggregate


@pytest.mark.order(99)
def test_chromnn_cpu(momics_path: str):
    mom = momics.Momics(momics_path)

    features = "bw2"
    features_size = 1024 + 1
    target = "bw3"
    stride = 48
    target_size = 24
    batch_size = 1000
    bins = mom.bins(width=features_size, stride=stride, cut_last_bin_out=True)
    bins2 = bins.copy()
    bins2.Start = bins2.Start + features_size // 2 - target_size // 2
    bins2.End = bins2.Start + target_size

    ### MANUAL TESTING, WITHOUT MOMICSDATASET ###
    X_train = query.MomicsQuery(mom, bins).query_tracks(tracks=[features]).coverage[features]
    X_train = np.array([X_train[chrom] for chrom in X_train.keys()])
    Y_train = query.MomicsQuery(mom, bins2).query_tracks(tracks=[target]).coverage[target]
    Y_train = np.array([Y_train[chrom] for chrom in Y_train.keys()])

    train_dataset = (
        tf.data.Dataset.from_tensor_slices(({features: X_train}, {target: Y_train}))
        .shuffle(buffer_size=len(X_train))
        .batch(batch_size)
        .prefetch(1)
        .repeat()
    )

    input = {features: layers.Input(shape=(features_size, 1), name=features)}
    output = {target: layers.Dense(target_size, activation="linear", name=target)}
    model = nn.ChromNN(input, output).model
    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.001), loss="mse")

    model.fit(
        train_dataset,
        epochs=10,
        steps_per_epoch=len(X_train) // batch_size,
    )

    bb = mom.bins(width=features_size, stride=5, cut_last_bin_out=True)["I", 0:50000]
    dat = query.MomicsQuery(mom, bb).query_tracks(tracks=["bw2"]).coverage["bw2"]
    dat = np.array(list(dat.values()))
    bb2 = bb.copy()
    bb2.Start = bb2.Start + features_size // 2 - target_size // 2
    bb2.End = bb2.Start + target_size
    chrom_sizes = {chrom: length for chrom, length in zip(mom.chroms().chrom, mom.chroms().length)}
    keys = [f"{chrom}:{start}-{end}" for chrom, start, end in zip(bb2.Chromosome, bb2.Start, bb2.End)]

    predictions = model.predict(dat)[target]
    res = {f"f{features_size}_s{stride}_t{target_size}": {k: None for k in keys}}
    for i, key in enumerate(keys):
        res[f"f{features_size}_s{stride}_t{target_size}"][key] = predictions[i]

    res = aggregate(res, bb2, chrom_sizes, type="mean", prefix="prediction")


def test_mae_cor():
    y_true = tf.constant([0.0, 1.0, 2.0, 3.0, 4.0], dtype=tf.float32)
    y_pred = tf.constant([0.0, 1.0, 2.0, 3.0, 4.0], dtype=tf.float32)

    nn.mae_cor(
        tf.reshape(y_true, [1, 5]),
        tf.reshape(y_pred, [1, 5]),
    )

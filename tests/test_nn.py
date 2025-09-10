import pytest

import numpy as np
import tensorflow as tf
from momics import momics
from momics import utils
from momics import query
from momics import nn
from momics import aggregate
from momics import attribution
from tensorflow.keras import layers  # type: ignore


## Deactivate GPU
tf.config.set_visible_devices([], "GPU")


@pytest.mark.order(99)
def test_chromnn_cpu():
    ## Initial vars
    mom = momics.Momics("tests_data/test.momics")
    features = "ATAC"
    features_size = 8192
    target = "SCC1"
    stride = 48
    target_size = 512
    batch_size = 100
    bins = mom.bins(width=features_size, stride=stride, cut_last_bin_out=True).sample(100)
    bins2 = bins.copy()
    bins2.Start = bins2.Start + features_size // 2 - target_size // 2
    bins2.End = bins2.Start + target_size

    ## Make tf reproducible
    tf.random.set_seed(42)
    np.random.seed(42)

    ## Manual dataset
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

    # Train model
    input = {features: layers.Input(shape=(features_size, 1), name=features)}
    output = {target: layers.Reshape((target_size,), name=target)}
    model = nn.ChromNN(input, output).model
    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.001), loss="mse")
    model.fit(train_dataset, epochs=2, steps_per_epoch=len(X_train) // batch_size)

    ## Predict
    bb = mom.bins(width=features_size, stride=5, cut_last_bin_out=True)["I", 0:100]
    dat = query.MomicsQuery(mom, bb).query_tracks(tracks=[features]).coverage[features]
    dat = np.array([dat[chrom] for chrom in dat.keys()])
    bb2 = bb.extend(-(features_size - target_size) // 2)
    chrom_sizes = {chrom: length for chrom, length in zip(mom.chroms().chrom, mom.chroms().length)}
    keys = [f"{chrom}:{start}-{end}" for chrom, start, end in zip(bb2.Chromosome, bb2.Start, bb2.End)]
    predictions = model.predict(dat)[target]
    res = {f"f{features_size}_s{stride}_t{target_size}": {k: None for k in keys}}
    for i, key in enumerate(keys):
        res[f"f{features_size}_s{stride}_t{target_size}"][key] = predictions[i]

    res = aggregate.aggregate(res, chrom_sizes, type="mean", prefix="prediction")
    assert len(res["f8192_s48_t512"]) == 17
    assert len(res["f8192_s48_t512"]["I"]) == 230218


@pytest.mark.order(99)
def test_chromnn_attribution():

    ## Initial vars
    mom = momics.Momics("tests_data/test.momics")
    features = "nucleotide"
    features_size = 128
    target = "ATAC"
    stride = 48
    target_size = 4
    batch_size = 100
    bins = mom.bins(width=features_size, stride=stride, cut_last_bin_out=True).sample(1000)
    bins2 = bins.copy()
    bins2.Start = bins2.Start + features_size // 2 - target_size // 2
    bins2.End = bins2.Start + target_size

    ## Make tf reproducible
    tf.random.set_seed(42)
    np.random.seed(42)

    ## Manual dataset
    seqs = list(query.MomicsQuery(mom, bins).query_sequence().seq[features].values())
    X_train = utils.one_hot_encode(seqs)
    Y_train = query.MomicsQuery(mom, bins2).query_tracks(tracks=[target]).coverage[target]
    Y_train = np.array([Y_train[chrom] for chrom in Y_train.keys()])
    train_dataset = (
        tf.data.Dataset.from_tensor_slices((X_train, Y_train))
        .shuffle(buffer_size=len(X_train))
        .batch(batch_size)
        .prefetch(1)
        .repeat()
    )
    seqmodel = nn.Basenji(input_size=features_size, output_size=target_size).model
    seqmodel.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.001), loss="mse")
    seqmodel.fit(train_dataset, epochs=2, steps_per_epoch=len(X_train) // batch_size)

    ## Test saliency
    vp_width = 100
    seq = X_train[0:1, :, :]
    sal = attribution.get_saliency(seqmodel, seq)
    assert len(sal) == 2
    assert sal[0].shape == (features_size,)
    assert sal[1].shape == (features_size, 4)

    ## Test ISM
    ism = attribution.get_ISM(seqmodel, seq, viewpoint_width=vp_width, batch=10)
    assert len(ism) == 3
    assert len(ism[0]) == vp_width
    assert ism[1].shape == (vp_width, 4)
    assert ism[2].shape == (vp_width,)

    ## Test full attribution
    attr = attribution.attribution(
        mom, seqmodel, track_name=target, chromosome="I", centerpoint=512, viewpoint_width=vp_width, batch=10
    )
    assert len(attr) == 8
    assert attr[0].shape == (1, features_size, 4)
    assert attr[1].shape == (vp_width,)
    assert attr[2].shape == (1, target_size)
    assert attr[3].shape == (features_size,)
    assert len(attr[4]) == vp_width
    assert attr[5].shape == (vp_width, 4)
    assert attr[6].shape == (vp_width,)
    assert attr[7].shape == (vp_width,)


def test_rc_augmentation():
    inputs = {"nucleotide": tf.expand_dims(tf.convert_to_tensor(utils.one_hot_encode("ATTGGGCC")), axis=0)}
    outputs = {"ATAC": tf.expand_dims(tf.convert_to_tensor(np.arange(8).astype("float32")), axis=0)}
    tf.random.set_seed(44)
    np.random.seed(44)
    augm = nn.tf_rc_augmentation(inputs, outputs)
    inp_rc = augm[0]["nucleotide"][0:1]
    out_rc = augm[1]["ATAC"][0:1]
    assert utils.one_hot_decode(inp_rc) == ["GGCCCAAT"]
    assert np.array_equal(out_rc.numpy(), np.array([[[7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0]]]))


def test_cor():
    y_true = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    y_pred = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    correlation = nn.cor(y_true, y_pred)
    assert np.isclose(correlation.numpy(), 1.0)

    y_pred = np.array([[3, 2, 1], [6, 5, 4]], dtype=np.float32)
    correlation = nn.cor(y_true, y_pred)
    assert np.isclose(correlation.numpy(), -1.0)


def test_loss_mae_cor():
    y_true = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    y_pred = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    loss = nn.loss_mae_cor(y_true, y_pred, alpha=0)  # loss is just correlation
    assert np.isclose(loss.numpy(), 0.0, atol=1e-6)
    loss = nn.loss_mae_cor(y_true, y_pred, alpha=1)  # loss is just mae
    assert np.isclose(loss.numpy(), 0.0, atol=1e-6)
    loss = nn.loss_mae_cor(y_true, y_pred, alpha=0.75)  # loss is both mae (75%) and cor (25%)
    assert np.isclose(loss.numpy(), 0.0, atol=1e-6)

    y_pred = np.array([[3, 2, 1], [6, 5, 4]], dtype=np.float32)
    loss = nn.loss_mae_cor(y_true, y_pred, alpha=0)  # loss is just correlation
    assert np.isclose(loss.numpy(), 2.0, atol=1e-6)
    loss = nn.loss_mae_cor(y_true, y_pred, alpha=1)  # loss is just mae
    assert np.isclose(loss.numpy(), 1.333333333, atol=1e-6)
    loss = nn.loss_mae_cor(y_true, y_pred, alpha=0.75)  # loss is both mae (75%) and cor (25%)
    assert np.isclose(loss.numpy(), 1.5, atol=1e-6)

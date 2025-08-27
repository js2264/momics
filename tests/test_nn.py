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


@pytest.mark.order(99)
def test_chromnn_cpu():

    ## Initial vars
    mom = momics.Momics("tests_data/yeast_CNN_data.momics")
    features = "ATAC_rescaled"
    features_size = 128
    target = "MNase_rescaled"
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
    output = {target: layers.Dense(target_size, activation="linear", name=target)}
    model = nn.ChromNN(input, output).model
    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.001), loss="mse")
    model.fit(train_dataset, epochs=2, steps_per_epoch=len(X_train) // batch_size)

    ## Predict
    bb = mom.bins(width=features_size, stride=5, cut_last_bin_out=True)["I", 0:50000]
    dat = query.MomicsQuery(mom, bb).query_tracks(tracks=[features]).coverage[features]
    dat = np.array([dat[chrom] for chrom in dat.keys()])
    bb2 = bb.copy()
    bb2.Start = bb2.Start + features_size // 2 - target_size // 2
    bb2.End = bb2.Start + target_size
    chrom_sizes = {chrom: length for chrom, length in zip(mom.chroms().chrom, mom.chroms().length)}
    keys = [f"{chrom}:{start}-{end}" for chrom, start, end in zip(bb2.Chromosome, bb2.Start, bb2.End)]
    predictions = model.predict(dat)[target]
    res = {f"f{features_size}_s{stride}_t{target_size}": {k: None for k in keys}}
    for i, key in enumerate(keys):
        res[f"f{features_size}_s{stride}_t{target_size}"][key] = predictions[i]

    res = aggregate.aggregate(res, chrom_sizes, type="mean", prefix="prediction")
    assert len(res["f128_s48_t4"]) == 17
    assert np.sum(res["f128_s48_t4"]["I"]) > 0


@pytest.mark.order(99)
def test_chromnn_attribution():

    ## Initial vars
    mom = momics.Momics("tests_data/yeast_CNN_data.momics")
    features = "nucleotide"
    features_size = 128
    target = "MNase_rescaled"
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
    assert attr[0].shape == (1, features_size, 5)
    assert attr[1].shape == (vp_width,)
    assert attr[2].shape == (1, target_size)
    assert attr[3].shape == (features_size,)
    assert len(attr[4]) == vp_width
    assert attr[5].shape == (vp_width, 4)
    assert attr[6].shape == (vp_width,)
    assert attr[7].shape == (vp_width,)

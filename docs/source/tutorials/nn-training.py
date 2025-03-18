import momics
from momics.momics import Momics
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers  # type: ignore
from pathlib import Path  # type: ignore
from tensorflow.keras.callbacks import CSVLogger, EarlyStopping, ReduceLROnPlateau  # type: ignore
from momics.aggregate import aggregate
from momics import utils as mutils
from momics.chromnn import ChromNN

# Fetch data from the momics repository
repo = Momics("yeast_CNN_data.momics")
features = "MNase_rescaled"
features_size = 1024 + 1
target = "ATAC_rescaled"
stride = 48
target_size = 24
batch_size = 1000
bins = repo.bins(width=features_size, stride=stride, cut_last_bin_out=True)
bins = bins.subset(lambda x: x.Chromosome != "XVI")


def one_hot_encode(seq) -> np.ndarray:
    seq = seq.upper()
    encoding_map = {"A": [1, 0, 0, 0], "T": [0, 1, 0, 0], "C": [0, 0, 1, 0], "G": [0, 0, 0, 1]}
    oha = np.zeros((len(seq), 4), dtype=int)
    for i, nucleotide in enumerate(seq):
        oha[i] = encoding_map[nucleotide]

    return oha


# Split data into training and validation sets
bins2 = bins.copy()
bins2.Start = bins2.Start + features_size // 2 - target_size // 2
bins2.End = bins2.Start + target_size
bins_split, bins_test = mutils.split_ranges(bins, 0.8, shuffle=True)
bins_train, bins_val = mutils.split_ranges(bins_split, 0.8, shuffle=True)
bins2_split, bins2_test = mutils.split_ranges(bins2, 0.8, shuffle=True)
bins2_train, bins2_val = mutils.split_ranges(bins2_split, 0.8, shuffle=True)
X_train = momics.query.MomicsQuery(repo, bins_train).query_tracks(tracks=[features]).coverage[features]
# X_train = momics.query.MomicsQuery(repo, bins_train).query_sequence().seq["nucleotide"]
X_train = np.array([X_train[chrom] for chrom in X_train.keys()])
# X_train = np.array([one_hot_encode(seq) for seq in X_train])
Y_train = momics.query.MomicsQuery(repo, bins2_train).query_tracks(tracks=[target]).coverage[target]
Y_train = np.array([Y_train[chrom] for chrom in Y_train.keys()])
X_val = momics.query.MomicsQuery(repo, bins_val).query_tracks(tracks=[features]).coverage[features]
# X_val = momics.query.MomicsQuery(repo, bins_val).query_sequence().seq["nucleotide"]
X_val = np.array([X_val[chrom] for chrom in X_val.keys()])
# X_val = np.array([one_hot_encode(seq) for seq in X_val])
Y_val = momics.query.MomicsQuery(repo, bins2_val).query_tracks(tracks=[target]).coverage[target]
Y_val = np.array([Y_val[chrom] for chrom in Y_val.keys()])
X_test = momics.query.MomicsQuery(repo, bins_test).query_tracks(tracks=[features]).coverage[features]
# X_test = momics.query.MomicsQuery(repo, bins_test).query_sequence().seq["nucleotide"]
X_test = np.array([X_test[chrom] for chrom in X_test.keys()])
# X_test = np.array([one_hot_encode(seq) for seq in X_test])
Y_test = momics.query.MomicsQuery(repo, bins2_test).query_tracks(tracks=[target]).coverage[target]
Y_test = np.array([Y_test[chrom] for chrom in Y_test.keys()])

train_dataset = (
    tf.data.Dataset.from_tensor_slices((X_train, Y_train))
    .shuffle(buffer_size=len(X_train))
    .batch(batch_size)
    .prefetch(1)
    .repeat()
)
val_dataset = tf.data.Dataset.from_tensor_slices((X_val, Y_val)).batch(batch_size).prefetch(1)

# ---------------------------------

# Define CNN
input = layers.Input(shape=(features_size, 4 if features == "seq" else 1))
output = layers.Dense(target_size, activation="linear")
model = ChromNN(input, output).model
model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.001), loss="mse")
model.summary()

# ---------------------------

# Model training

callbacks_list = [
    CSVLogger(Path(".chromnn", "epoch_data.csv")),
    EarlyStopping(monitor="val_loss", patience=8, min_delta=1e-5, restore_best_weights=True),
    ReduceLROnPlateau(monitor="val_loss", factor=0.1, patience=8 // 2, min_lr=0.1 * 0.001),
]

model.fit(
    train_dataset,
    validation_data=val_dataset,
    epochs=30,
    callbacks=callbacks_list,
    steps_per_epoch=len(X_train) // batch_size,
)

# ---------------------------

# Model evaluation
test_dataset = tf.data.Dataset.from_tensor_slices((X_test, Y_test)).batch(batch_size).prefetch(4)
model.evaluate(test_dataset)
model.save("simple_chromnn.keras")

# Reopen the model
# model = tf.keras.models.load_model("simple_chromnn.keras")
# model.summary()

# ---------------------------

# Model prediction
bb = repo.bins(width=features_size, stride=5, cut_last_bin_out=True)["XVI", 0:300000]
dat = momics.query.MomicsQuery(repo, bb).query_tracks(tracks=["MNase_rescaled"]).coverage["MNase_rescaled"]
# dat = momics.query.MomicsQuery(repo, bb).query_sequence().seq["nucleotide"]
dat = np.array(list(dat.values()))
# dat = np.array([one_hot_encode(seq) for seq in dat])

# Export predictions as a bigwig
bb2 = bb.copy()
bb2.Start = bb2.Start + features_size // 2 - target_size // 2
bb2.End = bb2.Start + target_size
chrom_sizes = {chrom: length for chrom, length in zip(repo.chroms().chrom, repo.chroms().length)}
keys = [f"{chrom}:{start}-{end}" for chrom, start, end in zip(bb2.Chromosome, bb2.Start, bb2.End)]

predictions = model.predict(dat)
res = {f"atac_f{features_size}_s{stride}_t{target_size}_finals5_2": {k: None for k in keys}}
for i, key in enumerate(keys):
    res[f"atac_f{features_size}_s{stride}_t{target_size}_finals5_2"][key] = predictions[i]

res = aggregate(res, bb2, chrom_sizes, type="mean", prefix="prediction")

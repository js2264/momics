import momics
from momics.momics import Momics
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers  # type: ignore
from pathlib import Path  # type: ignore
from tensorflow.keras.callbacks import CSVLogger, EarlyStopping, ReduceLROnPlateau  # type: ignore
from momics.aggregate import aggregate
from momics import utils as mutils
from momics.chromnn import Basenji

# Fetch data from the momics repository
repo = Momics("yeast_CNN_data.momics")
repo.tracks()
features = "seq"
features_size = 4096 + 1
target = "atac_rescaled"
stride = 256
target_size = 48
batch_size = 500
bins = repo.bins(width=features_size, stride=stride, cut_last_bin_out=True)
bins = bins.subset(lambda x: x.Chromosome != "XVI")

# Split data into training and validation sets
bins2 = bins.copy()
bins2.Start = bins2.Start + features_size // 2 - target_size // 2
bins2.End = bins2.Start + target_size
bins_split, bins_test = mutils.split_ranges(bins, 0.8, shuffle=False)
bins_train, bins_val = mutils.split_ranges(bins_split, 0.8, shuffle=False)
bins2_split, bins2_test = mutils.split_ranges(bins2, 0.8, shuffle=False)
bins2_train, bins2_val = mutils.split_ranges(bins2_split, 0.8, shuffle=False)

# Import data
X_train = momics.query.MomicsQuery(repo, bins_train).query_sequence().seq["nucleotide"]
X_train = np.array([X_train[chrom] for chrom in X_train.keys()])
X_train = np.array([mutils.one_hot_encode(seq) for seq in X_train])
Y_train = momics.query.MomicsQuery(repo, bins2_train).query_tracks(tracks=[target]).coverage[target]
Y_train = np.array([Y_train[chrom] for chrom in Y_train.keys()])
X_val = momics.query.MomicsQuery(repo, bins_val).query_sequence().seq["nucleotide"]
X_val = np.array([X_val[chrom] for chrom in X_val.keys()])
X_val = np.array([mutils.one_hot_encode(seq) for seq in X_val])
Y_val = momics.query.MomicsQuery(repo, bins2_val).query_tracks(tracks=[target]).coverage[target]
Y_val = np.array([Y_val[chrom] for chrom in Y_val.keys()])
X_test = momics.query.MomicsQuery(repo, bins_test).query_sequence().seq["nucleotide"]
X_test = np.array([X_test[chrom] for chrom in X_test.keys()])
X_test = np.array([mutils.one_hot_encode(seq) for seq in X_test])
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
model = Basenji(input, output).model
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

# Model prediction
bb = repo.bins(width=features_size, stride=12, cut_last_bin_out=True)["XVI", 0:300000]
dat = momics.query.MomicsQuery(repo, bb).query_sequence().seq["nucleotide"]
dat = np.array(list(dat.values()))
dat = np.array([mutils.one_hot_encode(seq) for seq in dat])

# Export predictions as a bigwig
bb2 = bb.copy()
bb2.Start = bb2.Start + features_size // 2 - target_size // 2
bb2.End = bb2.Start + target_size
chrom_sizes = {chrom: length for chrom, length in zip(repo.chroms().chrom, repo.chroms().length)}
keys = [f"{chrom}:{start}-{end}" for chrom, start, end in zip(bb2.Chromosome, bb2.Start, bb2.End)]

predictions = model.predict(dat)
res = {f"atac-from-seq_f{features_size}_s{stride}_t{target_size}": {k: None for k in keys}}
for i, key in enumerate(keys):
    res[f"atac-from-seq_f{features_size}_s{stride}_t{target_size}"][key] = predictions[i]

res = aggregate(res, bb2, chrom_sizes, type="mean", prefix="prediction")


#######################################################################
#######################################################################
#######################################################################
#######################################################################

##### USING MOMICS DATASET

import importlib
import momics
from momics.momics import Momics
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers  # type: ignore
from pathlib import Path  # type: ignore
from tensorflow.keras.callbacks import CSVLogger, EarlyStopping, ReduceLROnPlateau  # type: ignore
from momics.aggregate import aggregate
from momics import utils as mutils
from momics import chromnn
from momics.dataset import MomicsDataset

# Fetch data from the momics repository
repo = Momics("yeast_CNN_data.momics")
repo.tracks()
features = "nucleotide"
features_size = 2048 + 1
target = "mnase_rescaled"
stride = 48

bins = repo.bins(width=features_size, stride=stride, cut_last_bin_out=True)
bins = bins.subset(lambda x: x.Chromosome != "XVI")
bins_split, bins_test = mutils.split_ranges(bins, 0.8, shuffle=False)
bins_train, bins_val = mutils.split_ranges(bins_split, 0.8, shuffle=False)

target_size = 128
batch_size = 2000
train_dataset = (
    MomicsDataset(repo, bins_train, features, target, target_size=target_size, batch_size=batch_size)
    # .shuffle(10)
    .prefetch(2).repeat()
)
val_dataset = MomicsDataset(repo, bins_val, features, target, target_size=target_size, batch_size=batch_size).repeat()
test_dataset = MomicsDataset(repo, bins_test, features, target, target_size=target_size, batch_size=batch_size)

# ---------------------------------

# Define CNN
importlib.reload(chromnn)
input = layers.Input(shape=(features_size, 4 if features == "nucleotide" else 1))
output = layers.Dense(target_size, activation="linear")
model = chromnn.NucNN(input, output).model
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
    steps_per_epoch=len(bins_train) // batch_size,
    validation_steps=int(np.floor(len(bins_val) // batch_size)),
)

# ---------------------------

# Model prediction
bb = repo.bins(width=features_size, stride=12, cut_last_bin_out=True)["XVI", 0:300000]
dat = momics.query.MomicsQuery(repo, bb).query_sequence().seq["nucleotide"]
dat = np.array(list(dat.values()))
dat = np.array([mutils.one_hot_encode(seq) for seq in dat])
predictions = model.predict(dat)

# Export predictions as a bigwig
bb2 = bb.copy()
bb2.Start = bb2.Start + features_size // 2 - target_size // 2
bb2.End = bb2.Start + target_size
chrom_sizes = {chrom: length for chrom, length in zip(repo.chroms().chrom, repo.chroms().length)}
keys = [f"{chrom}:{start}-{end}" for chrom, start, end in zip(bb2.Chromosome, bb2.Start, bb2.End)]

res = {f"mnase-from-seq_f{features_size}_s{stride}_t{target_size}": {k: None for k in keys}}
for i, key in enumerate(keys):
    res[f"mnase-from-seq_f{features_size}_s{stride}_t{target_size}"][key] = predictions[i]

res = aggregate(res, bb2, chrom_sizes, type="mean", prefix="prediction")

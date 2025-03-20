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
from momics.dataset import MomicsDataset

# Fetch data from the momics repository
repo = Momics("yeast_CNN_data.momics")
features = "MNase_rescaled"
features_size = 4097 + 1
target = "ATAC_rescaled"
target_size = 24
stride = 48
batch_size = 2000
bins = repo.bins(width=features_size, stride=stride, cut_last_bin_out=True)
bins = bins.subset(lambda x: x.Chromosome != "XVI")
bins_split, bins_test = mutils.split_ranges(bins, 0.8, shuffle=True)
bins_train, bins_val = mutils.split_ranges(bins_split, 0.8, shuffle=True)

# ---------------------------------

# Create datasets
train_dataset = (
    MomicsDataset(repo, bins_train, features, target, target_size=target_size, batch_size=batch_size, silent=True)
    .shuffle(10)
    .prefetch(2)
    .repeat()
)
val_dataset = MomicsDataset(
    repo, bins_val, features, target, target_size=target_size, batch_size=batch_size, silent=True
).repeat()

# ---------------------------------

# Define CNN
input = layers.Input(shape=(features_size, 1))
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
    steps_per_epoch=len(bins_train) // batch_size,
    validation_steps=len(bins_val) // batch_size,
)

# ---------------------------

# Model evaluation
test_dataset = MomicsDataset(repo, bins_train, features, target, target_size=target_size, batch_size=batch_size)
model.evaluate(test_dataset)
model.save("simple_chromnn.keras")

# ---------------------------

# Model prediction
bb = repo.bins(width=features_size, stride=8, cut_last_bin_out=True)["XVI"]
dat = momics.query.MomicsQuery(repo, bb).query_tracks(tracks=["MNase_rescaled"]).coverage["MNase_rescaled"]
dat = np.array(list(dat.values()))

# # # Export predictions as a bigwig
bb2 = bb.copy()
bb2.Start = bb2.Start + features_size // 2 - target_size // 2
bb2.End = bb2.Start + target_size
chrom_sizes = {chrom: length for chrom, length in zip(repo.chroms().chrom, repo.chroms().length)}
keys = [f"{chrom}:{start}-{end}" for chrom, start, end in zip(bb2.Chromosome, bb2.Start, bb2.End)]

predictions = model.predict(dat)
res = {"atac3": {k: None for k in keys}}
for i, key in enumerate(keys):
    res["atac3"][key] = predictions[i]

res = aggregate(res, bb2, chrom_sizes, type="mean", prefix="prediction")

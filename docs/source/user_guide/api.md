# Introduction to `momics` API

```{danger}
This package is still under active development, and we make no promises
about the stability of any specific class, function, etc.
Pin versions if you're worried about breaking changes!
```

The `momics` package provides a Python API for creating and interacting with the `momics` data
model. This API is designed to be simple and intuitive, and is built on top of
the `tiledb` library.

## Creating `momics` repositories

The main entry point for the Python API is the `momics.Momics` class. This
class provides methods for creating and populating `momics` repositories.

```python
import momics
mom = momics.Momics("path_to_new_repo.momics")
```

The `momics.Momics` class has several attributes, including `path` and `cfg`.

```python
print(mom.path)
print(mom.cfg)
```

## Registering chromosomes

The first step after creating a `momics` repository is to register the chromosome
lengths. This can be done using the `ingest_chroms` method of the `momics.Momics` class.

```python
chr_lengths = {
    "I": 230218,
    "II": 813184,
    "III": 316620,
    "IV": 1531933,
    "V": 576874,
    "VI": 270161,
    "VII": 1090940,
    "VIII": 562643,
    "IX": 439888,
    "X": 745751,
    "XI": 666816,
    "XII": 1078177,
    "XIII": 924431,
    "XIV": 784333,
    "XV": 1091291,
    "XVI": 948066,
    "Mito": 85779
}
mom.ingest_chroms(chr_lengths, genome_reference = "S288c")
```

Once they are registered, chromosomes can be listed using `mom.chroms()` method.

```python
chroms = mom.chroms()
print(chroms)
```

## Populating `momics` repositories

Once the chromosome lengths have been registered, the `momics` repository can be
populated with `tracks`, `features` or `sequence` using the corresponding `ingest_*` method.

```python
# Ingest genome reference sequence
mom.ingest_sequence("path_to_genome.fa", threads = 18)

# Ingest genomic features
mom.ingest_features({
    "bed1": "path_to_bed1.bed",
    "bed2": "path_to_bed2.bed"
})

# Ingest genomic coverage tracks
mom.ingest_tracks({
    bw_a="path_to_bw_a.bw",
    bw_b="path_to_bw_b.bw",
    bw_c="path_to_bw_c.bw"
}, threads = 18)
```

The ingested data can be listed using the corresponding method.

```python
print(mom.sequence())
print(mom.features())
print(mom.tracks())
```

## Querying `momics` repositories

The `momics` package provides a dedicated `MomicsQuery` class,
to register query ranges, run queries and export results.

### Registering a query

```python
# Query by coordinates
from momics.query import MomicsQuery
import pyranges as pr

# Query over a single region
q = MomicsQuery(mom, "I:10-1000")

# Query using a BED file
gr = pr.read_bed("path_to_regions.bed")
q = MomicsQuery(mom, gr)
```

### Running a query

Once a query `q` is defined, it can be exectuted to extract data from
`sequence` and `tracks` tables.

```python
# Query sequences
q.query_sequence()
print(q.seq)

# Query specific tracks
q.query_tracks(tracks=["ATAC", "MNase"])
print(q.coverage)

# Query all available tracks
q.query_tracks()
print(q.to_df())
```

Both `query_*` methods provide a `threads` argument to parallelize the query
using the efficient tileDB storage backend.

```python
q.query_sequence(threads = 4)
q.query_tracks(threads = 4)
```

### Exporting query results

The query results can be coerced into generic bioinformatic data objects and
exported to output files using dedicated methods of the `MomicsQuery` class.

```python
# Coerce queried sequences as a SeqRecord object
q.to_SeqRecord()

# Export the queried scores as a json file
q.to_json("output.json")

# Export both sequences and scores as a npz file
q.to_npz("output.npz")

# Export as pandas DataFrame for further analysis
df = q.to_df()
df.to_csv("results.csv")
```

## Working with genomic bins and windows

`momics` provides utilities for creating genomic bins for systematic analysis:

```python
# Create genome-wide bins
bins = mom.bins(width=1000, stride=1000)  # Non-overlapping 1kb bins
bins_overlap = mom.bins(width=1000, stride=500)  # Overlapping bins

# Create bins for specific chromosomes
chr1_bins = mom.bins(width=1000, stride=1000)["I"]

# Create bins with specific properties
bins = mom.bins(
    width=2048,           # Window size
    stride=128,           # Step size
    cut_last_bin_out=True # Remove incomplete bins at chromosome ends
)

# Sample random bins for training
training_bins = bins.sample(10000)
```

## Data streaming and batch processing

For large-scale analysis, `momics` provides streaming interfaces:

```python
from momics.streamer import MomicsStreamer

# Create a data streamer
bins = mom.bins(width=1024, stride=128, cut_last_bin_out=True)
streamer = MomicsStreamer(
    mom,
    bins,
    features=["nucleotide", "atac"],
    batch_size=1000
)

# Each batch contains one-hot-encoded sequence and features for 1000 genomic windows
for batch in streamer:
    print(f"Processing batch {streamer.batch_index}/{streamer.num_batches}")
    nucleotide_data = batch["nucleotide"]  # Shape: (1000, 1024, 5)
    atac_data = batch["atac"]              # Shape: (1000, 1024, 1)

    # Process batch...
```

## Machine Learning integration

### Creating datasets for deep learning

```python
from momics.dataset import MomicsDataset

# Create a dataset for supervised learning
dataset = MomicsDataset(
    mom,
    bins,
    features=["nucleotide", "h3k27ac"],  # Input features
    target="rna_expression",             # Target variable
    target_size=128,                     # Output window size
    batch_size=32
)

# Create dataset with multiple targets
multi_target_dataset = MomicsDataset(
    mom,
    bins,
    features=["nucleotide"],
    target=["h3k27ac", "h3k4me3", "atac"],
    target_size=256,
    batch_size=16
)
```

### Using pre-built neural network architectures

```python
from momics import nn
from tensorflow.keras import layers  # type: ignore

# Use a ChromNN model for multi-modal input
inputs = {
    "nucleotide": layers.Input(shape=(1024, 5))
}
outputs = {
    "h3k27ac": layers.Dense(256, activation="linear"),
    "h3k4me3": layers.Dense(256, activation="linear"),
    "atac": layers.Dense(256, activation="linear")
}

model = nn.ChromNN(inputs, outputs).model
model.compile(
    optimizer="adam",
    loss={
        "h3k27ac": nn.loss_mae_cor,
        "h3k4me3": nn.loss_mae_cor,
        "atac": nn.loss_mae_cor
    },
    metrics={
        "h3k27ac": ["mae", nn.cor],
        "h3k4me3": ["mae", nn.cor],
        "atac": ["mae", nn.cor]
    })
model.fit(multi_target_dataset, epochs=10)
```

## Data management and repository operations

### Repository information and metadata

```python
# Get repository information
mom.path
mom.cfg

# List available data
mom.chroms()
mom.tracks()
mom.features()
mom.seq()
```

### Repository maintenance

```python
# Consolidate repository for optimal performance
mom.consolidate(vacuum=True)

# Create repository manifest
manifest = mom.manifest()

# Remove data
mom.remove_track("old_track")
```

### Copying and exporting data

```python
# Copy tracks to standard formats
mom.export_track("chip_seq", "output.bw")
mom.export_sequence("output.fa")
```

## Going further

- See the [full API reference](../api/index) for more information.
- Check out the `momics` [CLI quick guide](./cli) or the [full CLI reference](../cli/index) for more information.

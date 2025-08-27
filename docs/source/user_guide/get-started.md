# Get started with `momics`

```{danger}
This package is still under active development, and we make no promises
about the stability of any specific class, function, etc.
Pin versions if you're worried about breaking changes!
```

## Installation

With python `3.9` and higher, you can install `momics` from [PyPI](https://pypi.org/project/momics) using `pip`.

```shell
pip install momics
```

The dependencies will automatically be installed.

```{tip}
We highly recommend using the `conda` package manager to install scientific
packages like `momics`. To get `conda`, you can download either the
full [Anaconda](https://www.continuum.io/downloads) Python distribution
which comes with lots of data science software or the minimal
[Miniconda](http://conda.pydata.org/miniconda.html) distribution
which is just the standalone package manager plus Python.

In the latter case, you can install `momics` and all its dependencies as follows:

    conda install bioconda::momics

```

### Development installation

For the latest development version or to contribute to `momics`:

```shell
git clone https://github.com/js2264/momics.git
cd momics
pip install -e .
```

## Quick start

Here's a complete workflow to get you started with `momics`:

### 1. Create a new repository

```bash
# Create a new .momics repository from CLI
momics create my_experiment.momics
```

```py
# Or from Python
import momics
mom = momics.Momics("my_experiment.momics")
```

### 2. Set up chromosome information

```bash
# From command line using a chrom.sizes file
momics ingest chroms -g hg38 -f S288c.chrom.sizes my_experiment.momics
```

```py
# Or from Python
chr_lengths = {
    "chr1": 249250621,
    "chr2": 242193529,
    "chr3": 198295559,
    # ... other chromosomes
}
mom.ingest_chroms(chr_lengths, genome_reference="hg38")
```

### 3. Add genomic data

```bash
# Add reference genome sequence
momics ingest seq -f hg38.fa my_experiment.momics

# Add multiple signal tracks
momics ingest tracks \
    -f chip_h3k27ac=chip_h3k27ac.bw \
    -f atac_seq=atac_seq.bw \
    -f rna_seq=rna_seq.bw \
    my_experiment.momics

# Add genomic features/annotations
momics ingest features \
    -f genes=gencode.gtf \
    -f peaks=peaks.bed \
    -f enhancers=enhancers.bed \
    my_experiment.momics
```

```py
## Or from Python

# Ingest genome reference sequence
mom.ingest_sequence("hg38.fa", threads=18)

# Ingest genomic coverage tracks
mom.ingest_tracks({
    chip_h3k27ac="path_to_bw_a.bw",
    atac_seq="atac_seq.bw",
    rna_seq="rna_seq.bw"
}, threads=18)

# Ingest genomic features
mom.ingest_features({
    "genes": "genes.gtf",
    "peaks": "peaks.bed",
    "enhancers": "enhancers.bed"
})
```

### 4. Query your data

```bash
# Query a specific genomic region
momics query seq --coordinates "chr1:1000000-1001000" my_experiment.momics
momics query tracks --coordinates "chr1:1000000-1001000" my_experiment.momics

# Query using a BED file with multiple regions
momics query tracks --file regions_of_interest.bed my_experiment.momics

# Export results to files
momics query seq --file regions.bed -o sequences.fa my_experiment.momics
momics query tracks --file regions.bed -o results.json my_experiment.momics
```

```python
# Or from Python

import momics
from momics import MomicsQuery
import numpy as np
import pandas as pd

# Load your repository
mom = momics.Momics("my_experiment.momics")

# Query a genomic region
q = MomicsQuery(mom, "chr1:1000000-1002000")

# Get sequence data
q.query_sequence()
print("Sequence length:", len(q.seq["chr1"]))

# Get coverage tracks
q.query_tracks(tracks=["chip_h3k27ac", "atac_seq"])
coverage_data = q.coverage

# Convert to pandas DataFrame for analysis
df = q.to_df()
print(df.head())

# Calculate correlations between tracks
correlation = np.corrcoef(
    coverage_data["chip_h3k27ac"]["chr1"],
    coverage_data["atac_seq"]["chr1"]
)[0,1]
print(f"H3K27ac-ATAC correlation: {correlation:.3f}")
```

## Common workflows

### Analyzing chromatin landscape at enhancers

```python
# Query regions around transcription factor peaks
enhancers = mom.features("enhancers")
q = MomicsQuery(mom, enhancers)

# Get sequences and multiple histone marks
q.query_sequence()
q.query_tracks(tracks=["h3k27ac", "h3k4me1", "h3k4me3", "atac"])

# Export for downstream analysis
q.to_npz("tf_binding_data.npz")

# Calculate average profiles around peaks
import matplotlib.pyplot as plt

profiles = []
for region in q.coverage["h3k27ac"]:
    profiles.append(q.coverage["h3k27ac"][region])

mean_profile = np.mean(profiles, axis=0)
plt.plot(mean_profile)
plt.title("Average H3K27ac around enhancers")
plt.show()
```

### Preparing data for machine learning

```python
from momics.dataset import MomicsDataset

# Create genomic bins for training
bins = mom.bins(width=2048, stride=32, cut_last_bin_out=True)
training_bins = bins.sample(50000)

# Create dataset for sequence-to-epigenome prediction
dataset = MomicsDataset(
    mom,
    training_bins,
    features="nucleotide",          # DNA sequence input
    target="h3k27ac",               # Epigenomic targets
    target_size=128,
    batch_size=32
).map(lambda x, y: (x['nucleotide'], y['h3k27ac']))

# Train a neural network
from momics import nn
import tensorflow as tf

model = nn.Basenji(input_size=2048, output_size=128).model
model.compile(optimizer="adam", loss=nn.loss_mae_cor, metrics=[nn.cor])
model.fit(dataset, epochs=10, steps_per_epoch=1000)
```

### Genome-wide analysis

```python
# Create genome-wide bins
genome_bins = mom.bins(width=10000, stride=10000)  # 10kb non-overlapping bins

# Query all data genome-wide
q = MomicsQuery(mom, genome_bins)
q.query_tracks(tracks=["chip_h3k27ac", "atac_seq", "rna_seq"])

# Calculate genome-wide correlations
df = q.to_df()
correlation_matrix = df[["chip_h3k27ac", "atac_seq", "rna_seq"]].corr()
print(correlation_matrix)

# Find highly active regions
df["activity_score"] = (df["chip_h3k27ac"] + df["atac_seq"] + df["rna_seq"]) / 3
high_activity = df[df["activity_score"] > df["activity_score"].quantile(0.95)]

print(f"Found {len(high_activity)} highly active regions")
high_activity.to_bed("highly_active_regions.bed")
```

### Going further

- Check out the `momics` [API quick guide](api) or the [full API reference](../api/index) for more information.
- Check out the `momics` [CLI quick guide](./cli) or the [full CLI reference](../cli/index) for more information.
- Visit the [GitHub repository](https://github.com/your-org/momics) for issues and discussions
- Read more about TileDB data storage principles: [https://docs.tiledb.com/main](https://docs.tiledb.com/main)

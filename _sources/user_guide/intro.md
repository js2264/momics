# What is `momics`?

`momics` is a comprehensive Python package designed to simplify the storage, management, and analysis of multi-omics genomic data. Built on top of the high-performance [TileDB](https://tiledb.com/) storage engine, `momics` provides researchers with a unified framework for handling diverse genomic datasets efficiently.

## A file format for efficient storage of multi-omics data

At its core, `momics` introduces a standardized file format (`.momics`) that can store multiple types of genomic data in a single, compressed repository:

- **Genomic sequences**: Reference genomes and DNA sequences
- **Coverage tracks**: ChIP-seq, ATAC-seq, RNA-seq, MNase-seq, and other signal data
- **Genomic features**: Gene annotations, peaks, regulatory elements from BED/GFF files

The `.momics` repository format leverages `TileDB`'s columnar storage and compression algorithms to achieve:

- **Space efficiency**: Significantly reduced file sizes compared to traditional formats
- **Fast random access**: Quick retrieval of data from any genomic region
- **Scalability**: Handle datasets from single experiments to consortium-scale studies

```python
# Example: A single .momics file can contain everything
mom = momics.Momics("my_experiments.momics")

# Add reference genome sequence
mom.ingest_sequence("genome.fa")

# Add multiple signal tracks from different technologies
mom.ingest_tracks({
    "h3k27ac_tumor": "h3k27ac_tumor.bw",
    "h3k27ac_normal": "h3k27ac_normal.bw",
    "atac_tumor": "atac_tumor.bw",
    "atac_normal": "atac_normal.bw",
    "rna_tumor": "rna_tumor.bw",
    "rna_normal": "rna_normal.bw"
})

# Add genomic annotations
mom.ingest_features({
    "genes": "gencode.gtf",
    "enhancers": "enhancers.bed",
    "tumor_peaks": "tumor_specific_peaks.bed"
})
```

## An efficient query engine for multi-omics data

`momics` is designed to handle the complexity of modern genomic experiments. Its querying features include:

- A unified API for querying different genomic data modalities across large datasets;
- Batch querying and parallel processing for out-of-memory datasets;
- Built-in data validation and quality control checks to ensure data integrity;
- Integration with standard bioinformatics formats (FASTA, BED, JSON);
- NumPy arrays for machine learning workflows.

```python
# Example: Querying multi-omics data
from momics.query import MomicsQuery

# Query a specific genomic region
q = MomicsQuery(mom, "chr1:1000000-1002000")
q.query_tracks(tracks=["h3k27ac", "atac", "rna"])
q.query_sequence()

# Get results as a pandas DataFrame
df = q.to_df()

# Or export to files
q.to_json("region_data.json")
q.to_npz("region_data.npz")
```

## An integration layer for machine learning analysis of multi-omics data

`momics` bridges the gap between genomic data storage and modern machine learning workflows:

- Built-in data loaders and neural networks for integration with the Tensorflow framework;
- Support for both sequence and signal-based models;
- Custom loss functions for genomic prediction tasks;
- Attribution methods for model interpretability;

```python
# Example: Training a neural network on genomic data
from momics.dataset import MomicsDataset
from momics import nn

# Create a dataset to predict H3K27ac from sequence + ATAC-seq
dataset = MomicsDataset(
    mom,
    features=["nucleotide", "atac"],
    target="h3k27ac",
    feature_size=1024,
    target_size=128,
    batch_size=32
)

# Build and train a model
model = nn.ChromNN(
    inputs={"nucleotide": (1024, 5), "atac": (1024, 1)},
    outputs={"h3k27ac": 128}
).model

model.compile(optimizer="adam", loss=nn.loss_mae_cor)
model.fit(dataset, epochs=10)
```

## Key benefits

- **Simplified workflow**: One package for compressed and indexed storage, querying, and analysis of multi-omics data;
- **Performance**: Fast I/O and computation on large datasets;
- **Flexibility**: Support of standardized formats for I/O operations;
- **ML integration**: Seamless connection to ML frameworks.

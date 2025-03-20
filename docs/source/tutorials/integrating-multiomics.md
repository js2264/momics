---
title: Integrating multi-omics data with `momics`
toc-title: Table of contents
---

`momics` excels at managing multi-omics data. It allows you to store and
retrieve data from a variety of sources. Here, we will see how to
create, manage and query local repositories using `momics`.

## Creating a local `momics` repository

::::: {.cell execution_count="2"}
```python3
from momics.momics import Momics

## Creating repository
Momics("yeast_CNN_data.momics").remove()  # Purging existing files
repo = Momics("yeast_CNN_data.momics")  # Creating a brand new repository

## Registering chromosome sizes
## We will get chromosome sizes from a local fasta file.
from pyfaidx import Fasta

f = Fasta("/home/jaseriza/repos/momics/data/S288c.fa")
chrom_lengths = {chrom: len(seq) for chrom, seq in zip(f.keys(), f.values())}

repo.ingest_chroms(chrom_lengths, genome_version="S288c")
repo.chroms()

## Registering sequence
repo.ingest_sequence(f.filename)

## Registering genomic tracks
repo.ingest_tracks(
    {
        "ATAC": "/home/jaseriza/repos/momics/data/S288c_atac.bw",
        "MNase": "/home/jaseriza/repos/momics/data/S288c_mnase.bw",
    }
)
```

:::.cell-output
    momics :: INFO :: 2025-03-20 13:09:44,767 :: No cloud config found for momics.Consider populating `~/.momics.ini` file with configuration settings for cloud access.
    momics :: INFO :: 2025-03-20 13:09:44,780 :: Purged yeast_CNN_data.momics
    momics :: INFO :: 2025-03-20 13:09:44,784 :: Created yeast_CNN_data.momics
    momics :: INFO :: 2025-03-20 13:09:47,823 :: Genome sequence ingested in 3.0204s.
    momics :: INFO :: 2025-03-20 13:09:48,957 :: 2 tracks ingested in 1.1327s.
:::

:::.cell-output .cell-output-display execution_count="1"
    <momics.momics.Momics at 0x7724ecb3cf40>
:::
:::::

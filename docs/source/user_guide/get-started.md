# Get started with `momics`

## Installation

With python `3.8` and higher, you can install `momics`  from [PyPI](https://pypi.org/project/momics) using `pip`.

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

## Quick start

```bash
momics create my.momics
momics ingest chroms -f ~/genomes/S288c/S288c.chrom.sizes my.momics
momics ingest seq -f ~/genomes/S288c/S288c.fa my.momics
momics ingest tracks -f bw_a=mnase.bw -f bw_b=atac.bw -f bw_c=chip.bw my.momics
momics ingest features -f bed1=temp.bed my.momics
momics query seq --coordinates "I:10-1000" my.momics
momics query tracks --coordinates "I:10-1000" my.momics
```

## Going further

- Check out the `momics` [API quick guide](api) or the [full API reference](../api/index) for more information.
- Check out the `momics` [CLI quick guide](cli) or the [full CLI reference](../cli/index) for more information.
- Read more about TileDB data storage principles: [https://docs.tiledb.com/main](https://docs.tiledb.com/main)

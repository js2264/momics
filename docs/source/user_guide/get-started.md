# Get started

## Installation

With python `3.8` and higher, you can install `momics`  from [PyPI](https://pypi.org/project/momics) using `pip`.

```
pip install momics
```

The following requirements should be automatically installed:

- `tiledb`, `pyarrow`, `numpy`, `scipy`, `pandas`, `pyBigWig`.

```{tip}
We highly recommend using the `conda` package manager to install scientific 
packages like these. To get `conda`, you can download either the 
full [Anaconda](https://www.continuum.io/downloads) Python distribution 
which comes with lots of data science software or the minimal 
[Miniconda](http://conda.pydata.org/miniconda.html) distribution 
which is just the standalone package manager plus Python. 

In the latter case, you can install `momics` and all its dependencies as follows:

    conda install bioconda::momics

```

## Workflow

```
momics create hg19.momics 
momics add chroms hg19.chrom.sizes
momics add tracks a=sample1.bw b=sample2.bw c=sample3.bw hg19.momics 
momics tree hg19.momics
momics ls --table chroms hg19.momics
momics ls --table tracks hg19.momics
momics query --query "I:10-1000" hg19.momics
momics extract --out a.bw hg19.momics
```

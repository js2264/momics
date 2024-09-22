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

momics add chroms -f hg19.chrom.sizes hg19.momics
momics add seq -f hg19.fa hg19.momics 
momics add tracks -f bw_a=sample1.bw -f bw_b=sample2.bw -f bw_c=sample3.bw hg19.momics 

momics ls --table chroms hg19.momics
momics ls --table tracks hg19.momics

momics query seq --coordinates "I:10-1000" hg19.momics
momics query tracks --coordinates "I:10-1000" hg19.momics

momics export --track bw_b --output b_exported.bw hg19.momics
momics remove --track bw_c hg19.momics
```

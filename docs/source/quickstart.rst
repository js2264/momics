Quickstart
==========

Installation
------------

Install :py:mod:`momics`  from `PyPI <https://pypi.org/project/momics>`_ using :command:`pip`.

::

    $ pip install momics

Requirements:

- Python `3.8` and higher
- Python packages :py:mod:`tiledb`, :py:mod:`pyarrow`, :py:mod:`numpy`, :py:mod:`scipy`, :py:mod:`pandas`, :py:mod:`pyBigWig`.

We highly recommend using the `conda` package manager to install scientific 
packages like these. To get :command:`conda`, you can download either the 
full `Anaconda <https://www.continuum.io/downloads>`_ Python distribution 
which comes with lots of data science software or the minimal 
`Miniconda <http://conda.pydata.org/miniconda.html>`_ distribution 
which is just the standalone package manager plus Python. 

In the latter case, you can install :py:mod:`momics` and all its dependencies as follows:

::

    $ conda install bioconda::momics


Command line interface (CLI)
----------------------------

See the :ref:`CLI` reference for more information.


The :py:mod:`momics` package includes command-line tools for creating, querying and manipulating `.momics` files.

::

    $ momics create --name hg19.momics --chrom hg19.chrom.sizes
    $ momics ingest hg19.momics a=sample1.bw b=sample2.bw c=sample3.bw
    $ momics tree hg19.momics
    $ momics dump --table chroms hg19.momics
    $ momics dump --table tracks hg19.momics
    $ momics query --table a --query "I:10-1000" hg19.momics
    $ momics extract --table a --out a.bw hg19.momics


Python API
----------

See the :ref:`API` for more information.


Quickstart
==========

.. highlight:: console

Installation
------------

With python `3.8` and higher, you can install :py:mod:`momics`  from `PyPI <https://pypi.org/project/momics>`_ using :command:`pip`.

.. code-block:: bash

    pip install momics

The following requirements should be automatically installed:

- :py:mod:`tiledb`, :py:mod:`pyarrow`, :py:mod:`numpy`, :py:mod:`scipy`, :py:mod:`pandas`, :py:mod:`pyBigWig`.

.. We highly recommend using the `conda` package manager to install scientific 
.. packages like these. To get :command:`conda`, you can download either the 
.. full `Anaconda <https://www.continuum.io/downloads>`_ Python distribution 
.. which comes with lots of data science software or the minimal 
.. `Miniconda <http://conda.pydata.org/miniconda.html>`_ distribution 
.. which is just the standalone package manager plus Python. 

.. In the latter case, you can install :py:mod:`momics` and all its dependencies as follows:

.. .. code-block:: bash

..     conda install bioconda::momics


Command line interface (CLI)
----------------------------

The :py:mod:`momics` package includes command-line tools for creating, querying and manipulating `.momics` files.

.. code-block:: bash

    momics create hg19.momics 
    momics add chroms hg19.chrom.sizes
    momics add tracks a=sample1.bw b=sample2.bw c=sample3.bw hg19.momics 
    momics tree hg19.momics
    momics ls --table chroms hg19.momics
    momics ls --table tracks hg19.momics
    momics query --table a --query "I:10-1000" hg19.momics
    momics extract --table a --out a.bw hg19.momics

See the :ref:`CLI` reference for more information.

Python API
----------

See the :ref:`API` for more information.


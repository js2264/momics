![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fjs2264%2Fmomics%2Frefs%2Fheads%2Fdevel%2Fpyproject.toml)
[![PyPI - Version](https://img.shields.io/pypi/v/momics)](https://pypi.org/project/momics/)
[![Test & Doc](https://github.com/js2264/momics/actions/workflows/ci.yml/badge.svg)](https://github.com/js2264/momics/actions/workflows/ci.yml)
[![Codecov](https://img.shields.io/codecov/c/gh/js2264/momics)](https://app.codecov.io/gh/js2264/momics)
[![Black](https://img.shields.io/badge/style-black-black)](https://github.com/psf/black)

# momics

`momics` is both a revolutionary file format for efficient storage of ***multi-omics*** data, and a Python package designed to query and manipulate these files. The file format is specifically designed to handle genomic coverage tracks and sequences. The package provides an intuitive command-line interface (`CLI`) and a Python library `API` for bioinformatics workflows involving genomic data.

The `momics` package aims to facilitate:

* Creation: ingestion of genomic files into `momics` files;
* Query: sequential and range query patterns, with tabular and array retrieval;
* Scalability: cloud-native, out-of-core operations on the data;
* Distributability: data export in standard formats.

Follow `momics` development on [GitHub](https://github.com/js2264/momics).

```{toctree}
:caption: User Guide
:maxdepth: 1

user_guide/intro
user_guide/concepts
user_guide/get-started
user_guide/api
user_guide/cli
```

```{toctree}
:caption: Tutorials
:maxdepth: 1

tutorials/integrating-multiomics
tutorials/cloud-repos
tutorials/nn-training
```

```{toctree}
:caption: References
:maxdepth: 1

api/index
cli/index
```

```{toctree}
:caption: Changelog
:maxdepth: 1

changes
```

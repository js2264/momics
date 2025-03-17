![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fjs2264%2Fmomics%2Frefs%2Fheads%2Fdevel%2Fpyproject.toml)
[![PyPI - Version](https://img.shields.io/pypi/v/momics)](https://pypi.org/project/momics/)
[![Test & Doc](https://github.com/js2264/momics/actions/workflows/ci.yml/badge.svg)](https://github.com/js2264/momics/actions/workflows/ci.yml)
[![Codecov](https://img.shields.io/codecov/c/gh/js2264/momics)](https://app.codecov.io/gh/js2264/momics)
[![Black](https://img.shields.io/badge/style-black-black)](https://github.com/psf/black)

# momics

`momics` is both a revolutionary file format for efficient storage of ***multi-omics*** data, and a Python package designed to query and manipulate these files. The file format is specifically designed to handle genomic coverage tracks and sequences. The package provides an intuitive command-line interface (CLI) and a Python library for bioinformatics workflows involving genomic data.

## Install

You can install `momics` using `pip`:

```sh
pip install momics
```

Alternatively, clone this repository and install the package locally:

```sh
git clone https://github.com/js2264/momics.git
cd momics
pip install .
```

## Features

- Efficient genomic data storage: Store large genomic coverage tracks and genome reference sequences compactly.
- Multi-Range querying: Query multiple genomic regions simultaneously with high performance.
- Rich Python library: Directly access and manipulate genomic data using Python objects and methods.
- Full-fledged command-line interface (CLI): Perform common tasks such as adding new tracks, querying data, and extracting information directly from the shell.

## Usage

### CLI commands

- Add a track:

To ingest a `.bw` genomic coverage data into a momics repository, you can use the `ingest` command:

```sh
momics ingest tracks -f bw1=path/to/file.bw path/to/momics_repo
```

- Query genomic coverage:

You can query tracks using either UCSC-style coordinates or a BED file:

```sh
momics query tracks --coordinates "chr1:1-1000" path/to/momics_repo
momics query tracks --file path/to/file.bed path/to/momics_repo
```

### Python API

You can leverage advanced Python API to **load**, **query**, **stream** and use a momics repository in machine learning training pipelines.

```py
from momics.momics import Momics
from momics.streamer import MomicsStreamer
from momics.dataset import MomicsDataset

# Load a Momics repository
repo = Momics("path/to/momics_repo")

# Query tracks with coordinates
df = repo.query_tracks("<UCSC-style coordinates>")

# Stream tracks by batches of coordinates
coords = repo.bins(1024, stride=12, cut_last_bin_out=True)
stream = MomicsStreamer(repo, coords, batch_size=100)
for coverages in stream:
    #...

# Create a tensorflow training set from a momics repository
ds = (MomicsDataset(repo, coords, features="track_1", target="track_2", target_size=24, batch_size=100)
    .shuffle(10)
    .prefetch(1)
    .repeat()
)
ds.element_spec
# ((TensorSpec(shape=(None, 1024, 1), dtype=tf.float32, name=None),),
#  (TensorSpec(shape=(None, 24, 1), dtype=tf.float32, name=None),))
```

## Data format

`momics` uses a custom data format that combines genomic sequences and coverage tracks in a compressed and indexed form. The format allows for rapid access to any region of the genome and supports simultaneous querying of multiple genomic regions.

## Contributing

Contributions are welcome! Please submit pull requests or issues on the GitHub repository.

This project uses `black` to format code and `ruff` for linting. We also support `pre-commit` to ensure
these have been run. To configure your local environment, please install these development dependencies and set up
the commit hooks.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

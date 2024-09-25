# momics

`momics` is both a revolutionary file format for efficient storage of genomic data, and a Python package designed to query and manipulate these files. The file format is specifically designed to handle genomic coverage tracks and sequences. The package provides an intuitive command-line interface (CLI) and a Python library for bioinformatics workflows involving genomic data. 

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

### CLI Commands

- Add a track:

To ingest a `.bw` genomic coverage data into a momics repository, you can use the `add` command:

```sh
momics add tracks -f bw1=path/to/file.bw path/to/momics_repo
```

- Query genomic coverage: 

You can query tracks using either UCSC-style coordinates or a BED file:

```sh
momics query tracks --coordinates "chr1:1-1000" path/to/momics_repo
momics query tracks --file path/to/file.bed path/to/momics_repo
```

### Python API 

In Python, you can load and query a momics repository like this:

```py
from momics import Momics

# Load a Momics repository
repo = Momics("path/to/momics_repo")

# Query tracks with coordinates
df = repo.query_tracks("chr1:1-1000")
```

## Data Format

`momics` uses a custom data format that combines genomic sequences and coverage tracks in a compressed and indexed form. The format allows for rapid access to any region of the genome and supports simultaneous querying of multiple genomic regions.

## Contributing

Contributions are welcome! Please submit pull requests or issues on the GitHub repository.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

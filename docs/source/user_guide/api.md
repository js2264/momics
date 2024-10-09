# Python API

## Introduction

The `momics` package provides a Python API for creating and interacting with the `momics` data
model. This API is designed to be simple and intuitive, and is built on top of
the `tiledb` library.

## Creating `momics` repositories

The main entry point for the Python API is the `momics.Momics` class. This
class provides methods for creating and populating `momics` repositories.

```python
import momics
mom = momics.Momics("path_to_new_repo.momics")
```

The `momics.Momics` class has several attributes, including `path` and `cfg`.

```python
print(mom.path)
print(mom.cfg)
```

## Registering chromosomes

The first step after creating a `momics` repository is to register the chromosome
lengths. This can be done using the `add_chroms` method of the `momics.Momics` class.

```python
chr_lengths = {
    "I": 230218,
    "II": 813184,
    "III": 316620,
    "IV": 1531933,
    "V": 576874,
    "VI": 270161,
    "VII": 1090940,
    "VIII": 562643,
    "IX": 439888,
    "X": 745751,
    "XI": 666816,
    "XII": 1078177,
    "XIII": 924431,
    "XIV": 784333,
    "XV": 1091291,
    "XVI": 948066,
    "Mito": 85779
}
mom.add_chroms(chr_lengths, genome_reference = "S288c")
```

Once they are registered, chromosomes can be listed using `mom.chroms()` method.

```python
chroms = mom.chroms()
print(chroms)
```

## Populating `momics` repositories

Once the chromosome lengths have been registered, the `momics` repository can be
populated with `tracks`, `features` or `sequence` using the corresponding `add_*` method.

```python
# Ingest genome reference sequence
mom.add_sequence("path_to_genome.fa", threads = 18)

# Ingest genomic features
mom.add_features({
    "bed1": "path_to_bed1.bed",
    "bed2": "path_to_bed2.bed"
})

# Ingest genomic coverage tracks
mom.add_tracks({
    bw_a="path_to_bw_a.bw",
    bw_b="path_to_bw_b.bw",
    bw_c="path_to_bw_c.bw"
}, threads = 18)
```

The ingested data can be listed using the corresponding method.

```python
print(mom.sequence())
print(mom.features())
print(mom.tracks())
```

## Querying `momics` repositories

The `momics` package provides a dedicated `MultiRangeQuery` class,
to register query ranges, run queries and export results.

### Registering a query

```python
q = momics.MultiRangeQuery(mom, "I:10-1000")
```

### Running a query

Once a query `q` is defined, it can be exectuted to extract data from
`sequence` and `tracks` tables.

```python
q.query_sequence()
print(q.seq)

q.query_tracks()
print(q.coverage)
print(q.to_df())
```

Both `query_*` methods profide a `threads` argument to parallelize the query
using the efficient tileDB storage backend. By default, the number of threads
is set to all available threads.

```python
q.query_sequence(threads = 4)
q.query_tracks(threads = 4)
```

### Exporting query results

The query results can be coerced into generic bioinformatic data objects and
exported to output files using dedicated methods of the `MultiRangeQuery` class.

```python
# Coerce queried sequences as a SeqRecord object
q.to_SeqRecord()

# Export the queried scores as a json file
q.to_json("output.json")

# Export both sequences and scores as a npz file
q.to_npz("output.npz")
```

## Going further

- See the [full API reference](../api/index) for more information.
- Check out the `momics` CLI quick guide: [cli](cli)

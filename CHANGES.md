## Unreleased

## [v0.3.0](https://github.com/js2264/momics/releases/tag/0.4.0)

*Date: 2024-10-25*

### API changes

* `MultiRangeQuery` -> `MomicsQuery`.
* `manifest` method to output the whole configuration of a `Momics` instance.
* `consolidate` method to consolidate a `Momics` repository.
* No submodules are imported in `__init__.py` anymore.
* All export functions are now `Momics` methods.
* All ingestion methods for `Momics` now start with `ingest_*`.

### New features

* `aggregate` submodule to merge dictionaries of partial coverage tracks (over pyranges) into genome-wide tracks.
* Relatively basic `ChromNN` CNN `TensorFlow` model.
* `MomicsDataset` class to pass data to `TensorFlow` models.
* `MomicsStreamer` class to stream data from `Momics` repositories.
* `Momics` can ingest `pyranges` objects to store genomic features.
* Sequences, features and tracks can be extracted from a `Momics` repository.

### Enhancements

* CLI is now partially based on `cloup` to improve user experience.
* `seq(label = "...)` now returns the sequence for an entire chromosome.
* `tracks(label = "...)` now returns the genome-wide track.
* `MomicsQuery` queries now rely on `pyranges` for range queries.
* `MomicsQuery` queries can extract only a subset of the tracks.
* `add_track` and `remove` methods for `Momics` class.
* `to_npz` and `to_json` methods for `MomicsQuery` class.

### Maintenance

* Support jupyter notebooks in documentation.
* Add changelog.
* Improve docs.
* Logging system updates.
* CLI updates.
* Add `codecov` support.

### Bug fixes

* Ensure that queries are done per chromosome, even if provided ranges are stranded.
* All ranges are now 0-based half-open intervals, as in `pyranges` and BED files.
* Parallelization for queries only relies on `TileDB` internal system.
* Filters used in `Momics` tables.
* Removal of `Azure`-hosted repositories.

## [v0.3.0](https://github.com/js2264/momics/releases/tag/0.3.0)

*Date: 2024-09-26*

### New features

* `MultiRangeQuery` relies on `pybedtools` for range queries.
* Queries support both coordinates and multi-loci queries.
* Added `MultiRangeQuery` class for queries.
* New `bins` method for `Momics`.
* Repositories can ingest fasta files.
* Added `export` submodule.
* Added `MomicsConfig` class to manage cloud configuration settings.
* Added `Click`-based CLI commands.

### Enhancements

* Improve query filters.
* CLI commands have mutli-threading options.

### Maintenance

* Added `ruff` support.
* Changed doc theme.
* Added support for autoapi in documentation.
* Improved tests coverage.

### Bug fixes

* Fixed broken use of single-file temporary store in `create_from_unordered`.
* Added heuristic in pairix cload to prevent excessively large chunks. #92
* Added extra checks in `cload pairix` and `cload tabix`. #62, #75

## [v0.2.0](https://github.com/js2264/momics/releases/tag/0.2.0)

*Date: 2024-07-31*

### Enhancements

* Implement `Momics` class.

## [v0.1.0](https://github.com/js2264/momics/releases/tag/0.1.0)

*Date: 2024-07-29*

### Enhancements

* Initial prototype.

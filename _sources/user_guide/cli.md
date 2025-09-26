# Introduction to `momics` CLI

The `momics` package includes command-line tools for creating, querying and manipulating `.momics` files.


## Basic CLI usage

```shell
# Check version and get help
momics -v
momics --help

# Initiate a momics repository
momics create testCLI.momics

# Register chromosome lengths
momics ingest chroms -f S288c.chrom.sizes testCLI.momics

# Ingest genome reference sequence
momics ingest seq -f S288c.fa testCLI.momics

# Ingest multiple tracks with different names
momics ingest tracks -f bw_a=track1.bw -f bw_b=track2.bw -f bw_c=track3.bw testCLI.momics

# Batch track ingest from a directory
momics ingest bulk /path/to/bigwig/files.bw my_experiment.momics

# Ingest features from different file formats
momics ingest features -f bed1=regions1.bed -f bed2=regions2.bed testCLI.momics

# Print all created tables and arrays
momics tree testCLI.momics

# Generate a manifest of the repository configuration and timestamps
momics manifest -o manifest.json testCLI.momics

# Consolidate the repository to optimize storage and performance
momics consolidate --vacuum testCLI.momics

# Summary of each table
momics ls --table chroms testCLI.momics
momics ls --table tracks testCLI.momics
momics ls --table features testCLI.momics

# Perform queries
momics query seq --coordinates "I:10-1000" testCLI.momics
momics query seq --file regions1.bed -o out.fa testCLI.momics
momics query tracks --coordinates "I:10-1000" testCLI.momics
momics query tracks --file regions1.bed testCLI.momics
```

## Repository management

```shell
# Get detailed repository information
momics info my_experiment.momics
momics tree my_experiment.momics

# List available data with metadata
momics ls --table chroms my_experiment.momics
momics ls --table tracks my_experiment.momics

# Consolidate the repository to optimize storage and performance
momics consolidate --vacuum my_experiment.momics

# Delete a repository
momics delete my_experiment.momics
```

## Data export and format conversion

```shell
# Export individual tracks
momics cp --type track --label chip_h3k27ac --output h3k27ac.bw my_experiment.momics
momics cp --type track --label atac_seq --output atac.bw my_experiment.momics

# Export features in different formats
momics cp --type features --label genes --output genes.gtf my_experiment.momics
momics cp --type features --label peaks --output peaks.bed my_experiment.momics

# Export sequences
momics cp --type sequence --output genome.fa my_experiment.momics
momics cp --type sequence --coordinates "chr1:1000-2000" --output region.fa my_experiment.momics
```

## Going further

- See the [full CLI reference](../cli/index) for more information.

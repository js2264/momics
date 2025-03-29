import collections
from pathlib import Path
from typing import List, Tuple, Union

import numpy as np
import pyranges as pr
import pyBigWig
import pyfaidx


def _repo_exists(path, cfg) -> bool:
    x = cfg.vfs.is_dir(path)
    if x:
        return True
    else:
        return False


def _check_fasta_lengths(fasta, chroms) -> None:
    reference_lengths = dict(zip(chroms["chrom"], chroms["length"]))
    if isinstance(fasta, Path):
        fasta = fasta.name
    with pyfaidx.Fasta(fasta) as fa:
        lengths = {name: len(seq) for name, seq in fa.items()}
    if lengths != reference_lengths:
        raise Exception(f"{fa} file do not have identical chromomosome lengths.")


def _check_chr_lengths(bw_files, chroms) -> None:
    reference_lengths = dict(zip(chroms["chrom"], chroms["length"]))
    for file in list(bw_files.values()):
        if isinstance(file, Path):
            file = file.name
        with pyBigWig.open(file) as bw:
            lengths = bw.chroms()
            if lengths != reference_lengths:
                raise Exception(f"{file} files do not have identical chromomosome lengths.")


def _check_track_names(bw_files, tracks) -> None:
    labels = set(tracks["label"])
    for element in list(bw_files.keys()):
        if element in labels:
            raise ValueError(f"Provided label '{element}' already present in `tracks` table")


def _check_feature_names(features, sets) -> None:
    labels = set(sets["label"])
    for element in list(features.keys()):
        if element in labels:
            raise ValueError(f"Provided label '{element}' already present in `features` table")


def _check_feature_name(feature, features) -> None:
    labels = set(features["label"])
    if feature not in labels:
        raise ValueError(f"Provided feature name '{feature}' does not exist in `features` table")


def _check_track_name(track, tracks) -> None:
    labels = set(tracks["label"])
    if track not in labels:
        raise ValueError(f"Provided track name '{track}' does not exist in `tracks` table")


def get_chr_lengths(bw: Union[Path, str]) -> dict:
    """
    A simple wrapper around pyBigWig to get chromosome lengths from a bigwig file.

    Args:
        bw (Path): path to a bigwig file

    Returns:
        dict: Dictionary of chromosome lengths
    """
    if isinstance(bw, Path):
        bw = bw.name
    with pyBigWig.open(bw) as b:
        a = b.chroms()
    b.close()
    return a


def dict_to_bigwig(bw_dict: dict, output: Union[Path, str]) -> Path:
    """
    Write a dictionary of coverages to a bigwig file.
    The dictionary should have chromosome names as keys and per-base coverage as values.

    Args:
        bw_dict (dict): Dictionary of chromosome coverages
        output (Path): Path to output bigwig file

    Returns:
        Path to the output bigwig file

    Examples:
        >>> bw_dict = {'chr1': np.random.rand(1000), 'chr2': np.random.rand(2000)}
        >>> dict_to_bigwig(bw_dict, 'output.bw')
    """
    if isinstance(output, Path):
        output = output.name

    bw = pyBigWig.open(output, "w")
    header = [(chrom, len(coverage)) for chrom, coverage in bw_dict.items()]
    bw.addHeader(header)
    for chrom, coverage in bw_dict.items():
        values0 = np.float32(coverage)
        bw.addEntries(chrom, 0, values=values0, span=1, step=1)
    bw.close()

    return Path(output)


def parse_ucsc_coordinates(coords: Union[List, str]) -> pr.PyRanges:
    """
    Parse UCSC-style coordinates as a pr.PyRanges object. The coordinates should be in the format "chrom:start-end".

    Args:
        coords (str): A UCSC-style set of coordinates (e.g., "I:11-100").

    Returns:
        pr.PyRanges: A pr.PyRanges object.
    """
    if isinstance(coords, str):
        coords = [coords]

    coords_dict = collections.defaultdict(list)
    for coord in coords:
        try:
            chr_part, range_part = coord.split(":")
            start, end = range_part.split("-")
            start = int(start)
            end = int(end)
            coords_dict["chr"].append(chr_part)
            coords_dict["start"].append(start)
            coords_dict["end"].append(end)

        except ValueError as e:
            raise ValueError(f"Invalid start/end values in coordinate '{coord}'. " + "Start and end must be integers.") from e
        except Exception as e:
            raise ValueError(
                f"Invalid format for UCSC-style coordinate '{coord}'. " + "Expected format: 'chrom:start-end'."
            ) from e

    return pr.PyRanges(chromosomes=coords_dict["chr"], starts=coords_dict["start"], ends=coords_dict["end"])


def split_ranges(pyranges, ratio=0.8, shuffle=True) -> Tuple[pr.PyRanges, pr.PyRanges]:
    """
    Split a PyRanges object into two PyRanges objects based on a ratio.
    The first PyRanges object will contain the first `ratio` proportion of the
    ranges, and the second PyRanges object will contain the remaining ranges.

    Args:
        pyranges (pr.PyRanges): A PyRanges object.
        ratio (float): A float between 0 and 1.

    Returns:
        Tuple[pr.PyRanges, pr.PyRanges]: A tuple of two PyRanges objects.
    """
    df = pyranges.df
    if shuffle:
        df = pyranges.df.sample(frac=1, random_state=42).reset_index(drop=True)
    split_idx = int(len(df) * ratio)
    train_df = df.iloc[:split_idx]
    test_df = df.iloc[split_idx:]
    train_pyranges = pr.PyRanges(train_df)
    test_pyranges = pr.PyRanges(test_df)
    return train_pyranges, test_pyranges


def pyranges_to_bw(pyranges: pr.PyRanges, scores: np.ndarray, output: str) -> None:
    """
    Write a PyRanges object and corresponding scores to a BigWig file.
    The PyRanges object must have the same length as the first dimension of the scores array.
    The PyRanges object must have ranges of the same width as the second dimension of the scores array.

    Args:
        pyranges (pr.PyRanges): A PyRanges object.
        scores (np.ndarray): A 2D NumPy array of scores.
        output (str): Path to the output BigWig file.

    Returns:
        None
    """
    # Abort if output file already exists
    if Path(output).exists():
        raise FileExistsError(f"Output file '{output}' already exists")

    # Check that pyranges length is the same as scores dim 0
    if len(pyranges) != scores.shape[0]:
        raise ValueError("Length of PyRanges object must be the same as scores dimension 0")

    # Check that all pyranges widths are equal to the scores dim 1
    widths = pyranges.End - pyranges.Start
    if len(set(widths)) != 1:
        raise ValueError("All ranges must have the same width")
    if next(iter(widths)) != scores.shape[1]:
        raise ValueError("All ranges must have the same width as the second dimension of scores")

    # Save chrom sizes in header
    bw = pyBigWig.open(output, "w")
    # if there is only one chromosome, get its size and add it to the header
    if len(pyranges.Chromosome.unique()) == 1:
        chrom_size = pyranges.df["End"].max()
        bw.addHeader([(next(iter(pyranges.Chromosome)), chrom_size)])
    else:
        chrom_sizes = pyranges.df.groupby("Chromosome", observed=False)["End"].max().to_dict()
        chroms = list(chrom_sizes.keys())
        sizes = list(chrom_sizes.values())
        bw.addHeader(list(zip(chroms, sizes)))

    # Iterate over the PyRanges and write corresponding scores
    df = pyranges.df
    for i, (chrom, start, end) in enumerate(zip(df.Chromosome, df.Start, df.End)):
        score = scores[i]
        positions = list(range(start, end))
        bw.addEntries([chrom] * len(positions), positions, ends=[p + 1 for p in positions], values=score)

    # Step 4: Close the BigWig file
    bw.close()


def one_hot_encode(sequences, handle_non_standard=False, dtype=np.int8) -> np.ndarray:
    """
    Efficiently one-hot encode DNA sequences.

    Args:
        sequences (Union[str, List[str]]): A single DNA sequence or list of DNA sequences
        handle_non_standard (bool): If True, non-standard nucleotides (not A,T,G,C) will be
                                   encoded as [0,0,0,0]. If False, will raise a KeyError.
        dtype: NumPy data type for the output array (default: np.int8 to save memory)

    Returns:
        np.ndarray: A one-hot encoded array of shape (len(sequences), seq_length, 4) for multiple
                  sequences or (seq_length, 4) for a single sequence
    """
    # Handle single sequence
    single_input = isinstance(sequences, str)
    if single_input:
        sequences = [sequences]

    # Convert all sequences to uppercase
    sequences = [seq.upper() for seq in sequences]

    # Get max sequence length for the output array
    max_len = max(len(seq) for seq in sequences)

    # Initialize output array
    output = np.zeros((len(sequences), max_len, 4), dtype=dtype)

    # Define the mapping for standard nucleotides
    mapping = {"A": 0, "T": 1, "C": 2, "G": 3}

    # Process each sequence using NumPy vectorization
    for i, seq in enumerate(sequences):
        seq_len = len(seq)

        if handle_non_standard:
            # Process each nucleotide, skipping non-standard ones
            for j, nt in enumerate(seq):
                if nt in mapping:
                    output[i, j, mapping[nt]] = 1
        else:
            # Process only standard nucleotides, faster but can raise KeyError
            indices = np.array([mapping[nt] for nt in seq])
            output[i, np.arange(seq_len), indices] = 1

    # Return single sequence format if input was a single sequence
    if single_input:
        return output[0, : len(sequences[0])]

    return output

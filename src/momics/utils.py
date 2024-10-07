from pathlib import Path

import pandas as pd
import pybedtools
import pyBigWig
import pyfaidx


def get_chr_lengths(bw: Path) -> dict:
    """Parse bigwig header to extract chromosome lengths

    Args:
        bw (Path): path to a bigwig file

    Returns:
        dict: Dictionary of chromosome lengths
    """
    with pyBigWig.open(bw) as bw:
        a = bw.chroms()
    bw.close()
    return a


def _dict_to_bigwig(bw_dict: dict, output: Path):
    """Write a dictionary of coverages to a bigwig file

    Args:
        bw_dict (dict): Dictionary of chromosome coverages
        output (Path): Path to output bigwig file
    """
    bw = pyBigWig.open(output, "w")
    header = [(chrom, len(coverage)) for chrom, coverage in bw_dict.items()]
    bw.addHeader(header)
    for chrom, coverage in bw_dict.items():
        bw.addEntries(chrom, 1, span=1, step=1, values=coverage)
    bw.close()


def _check_fasta_lengths(fasta, chroms):
    reference_lengths = dict(zip(chroms["chrom"], chroms["length"]))
    with pyfaidx.Fasta(fasta) as fa:
        lengths = {name: len(seq) for name, seq in fa.items()}
    if lengths != reference_lengths:
        raise Exception(f"{fa} file do not have identical chromomosome lengths.")


def _check_chr_lengths(bw_files, chroms):
    reference_lengths = dict(zip(chroms["chrom"], chroms["length"]))
    for file in list(bw_files.values()):
        with pyBigWig.open(file) as bw:
            lengths = bw.chroms()
            if lengths != reference_lengths:
                raise Exception(f"{file} files do not have identical chromomosome lengths.")


def _check_track_names(bw_files, tracks):
    labels = set(tracks["label"])
    for element in list(bw_files.keys()):
        if element in labels:
            raise ValueError(f"Provided label '{element}' already present in `tracks` table")


def _check_feature_names(features, sets):
    labels = set(sets["label"])
    for element in list(features.keys()):
        if element in labels:
            raise ValueError(f"Provided label '{element}' already present in `features` table")


def _check_feature_name(feature, features):
    labels = set(features["label"])
    if feature not in labels:
        raise ValueError(f"Provided feature name '{feature}' does not exist in `features` table")


def _check_track_name(track, tracks):
    labels = set(tracks["label"])
    if track not in labels:
        raise ValueError(f"Provided track name '{track}' does not exist in `tracks` table")


def parse_ucsc_coordinates(coords: str) -> pd.DataFrame:
    """Parse UCSC-style coordinates as a pybedtools.BedTool object.

    Args:
        coords (str): A UCSC-style set of coordinates (e.g., "I:11-100"). Note
        that the coordinates are 1-based.

    Returns:
        pd.DataFrame: A pybedtools.BedTool object.
    """
    if isinstance(coords, str):
        coords = [coords]

    coord_strings = []
    for coord in coords:
        try:
            chr_part, range_part = coord.split(":")
            start, end = range_part.split("-")
            start = int(start)
            end = int(end)
            coord_strings.append(f"{chr_part} {start} {end}")

        except ValueError as e:
            raise ValueError(f"Invalid start/end values in coordinate '{coord}'. " + "Start and end must be integers.") from e
        except Exception as e:
            raise ValueError(
                f"Invalid format for UCSC-style coordinate '{coord}'. " + "Expected format: 'chrom:start-end'."
            ) from e

    return pybedtools.BedTool("\n".join(coord_strings), from_string=True)

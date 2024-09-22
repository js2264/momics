from pathlib import Path

import pandas as pd
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


def _check_fasta_lengths(fasta, chroms):
    reference_lengths = dict(zip(chroms["chr"], chroms["length"]))
    with pyfaidx.Fasta(fasta) as fa:
        lengths = {name: len(seq) for name, seq in fa.items()}
    if lengths != reference_lengths:
        raise Exception(f"{fa} file do not have identical chromomosome lengths.")


def _check_chr_lengths(bw_files, chroms):
    reference_lengths = dict(zip(chroms["chr"], chroms["length"]))
    for file in list(bw_files.values()):
        with pyBigWig.open(file) as bw:
            lengths = bw.chroms()
            if lengths != reference_lengths:
                raise Exception(
                    f"{file} files do not have identical chromomosome lengths."
                )


def _check_track_names(bw_files, tracks):
    labels = set(tracks["label"])
    for element in list(bw_files.keys()):
        if element in labels:
            raise ValueError(
                f"Provided label '{element}' already present in `tracks` table"
            )


def _check_chr_name(chr, chroms):
    if chr not in chroms["chr"].values:
        raise ValueError(f"{chr} chromosome not listed in `chroms` table.")


def _check_track_name(track, tracks):
    labels = set(tracks["label"])
    if track not in labels:
        raise ValueError(
            f"Provided track name '{track}' does not exist in `tracks` table"
        )


def import_bed_file(file_path: Path) -> pd.DataFrame:
    """Import a local bed file as a DataFrame

    Args:
        file_path (Path): Path to a local bed file

    Returns:
        pd.DataFrame: A pd.DataFrame with at least three columns (`chr`, `start`, `end`),
        and any other columns stored in the local bed file.
    """
    # Convert file_path to a Path object
    file_path = Path(file_path)

    try:
        # Read the file into a DataFrame
        df = pd.read_csv(file_path, sep="\t", header=None)

        # Check if the DataFrame has at least 3 columns
        if df.shape[1] < 3:
            print(f"The file {file_path} does not have at least 3 columns.")
            return False

        # Ensure that columns 2 and 3 are numeric (start and end coordinates)
        df[1] = pd.to_numeric(
            df[1], errors="coerce"
        )  # Convert start column to numeric, set errors to NaN
        df[2] = pd.to_numeric(
            df[2], errors="coerce"
        )  # Convert end column to numeric, set errors to NaN

        if df[1].isna().any() or df[2].isna().any():
            print(
                f"Columns 2 and 3 in {file_path} must contain valid numeric "
                + "values for start and end coordinates."
            )
            return False

    except Exception as e:
        raise e

    df.columns = ["chr", "start", "end"] + [col for col in df.columns[3:]]
    return df


def parse_ucsc_coordinates(coords: str) -> pd.DataFrame:
    """Parse UCSC-style coordinates as a DataFrame

    Args:
        coords (str): A UCSC-style set of coordinates (e.g., "I:11-100"). Note
        that the coordinates are 1-based.

    Returns:
        pd.DataFrame: A pd.DataFrame with at least three columns (`chr`, `start`, `end`),
        and any other columns stored in the local bed file.
    """
    if isinstance(coords, str):
        coords = [coords]

    chromosomes = []
    starts = []
    ends = []
    for coord in coords:
        try:
            chr_part, range_part = coord.split(":")
            start, end = range_part.split("-")
            start = int(start)
            end = int(end)
            chromosomes.append(chr_part)
            starts.append(start)
            ends.append(end)

        except ValueError:
            return (
                f"Error: Invalid start/end values in coordinate '{coord}'. "
                + "Start and end must be integers."
            )
        except Exception:
            return (
                f"Error: Invalid format for UCSC-style coordinate '{coord}'. "
                + "Expected format: 'chr:start-end'."
            )

    df = pd.DataFrame({"chr": chromosomes, "start": starts, "end": ends})

    return df

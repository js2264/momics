from typing import Optional, Tuple  

import numpy as np
import pyranges as pr
import logging
import tensorflow as tf

from .momics import Momics
from .multirangequery import MultiRangeQuery


class RangeDataLoader(tf.keras.utils.Sequence):
    """
    This class is implemented to train deep learning models, where the
    input data is a track or a sequence and the label is another track.
    The generator will iterate over the ranges in batches and extract
    the data and label for each range.

    Attributes
    ----------
    momics (Momics): a local `.momics` repository.
    ranges (dict): pr.PyRanges object.
    data (str): the name of the track to use as data
    label (str): the name of the track to use as label
    label_size (int): To which width should the label be centered
    """

    def __init__(
        self,
        momics: Momics,
        ranges: pr.PyRanges,
        batch_size: Optional[int],
        data: str,
        label: str,
        label_size: Optional[int] = None,
        silent: bool = False,
    ) -> None:
        """Initialize the RangeGenerator object.

        Args:
            momics (Momics): a Momics object
            ranges (pr.PyRanges): pr.PyRanges object
            batch_size (int): the batch size
            data (str): the name of the track to use as data
            label (str): the name of the track to use as label
            label_size (int): To which width should the label be centered
        """

        # Check that all ranges have the same width
        df = ranges.df
        widths = df.End - df.Start + 1
        if len(set(widths)) != 1:
            raise ValueError("All ranges must have the same width")

        self.momics = momics
        self.ranges = ranges
        if batch_size is None:
            batch_size = len(ranges)
        self.start = 0
        self.stop = len(ranges)
        self.batch_size = batch_size
        self.current = self.start
        self.silent = silent

        tr = momics.tracks()
        if data == "nucleotide":
            _ = momics.seq()

        if data not in list(tr["label"]) and data != "nucleotide":
            raise ValueError(f"Track {data} not found in momics repository.")
        if label not in list(tr["label"]):
            raise ValueError(f"Track {label} not found in momics repository.")

        self.data = data
        self.label = label

        if label_size is not None and label_size >= int(widths[0]):
            raise ValueError("Label center must be smaller than the range width.")
        self.label_size = label_size

    def __len__(self) -> int:
        return int(np.ceil(len(self.ranges) / self.batch_size))

    def __getitem__(self, idx) -> Tuple[np.ndarray, np.ndarray]:
        subrg = pr.PyRanges(self.ranges.df[idx * self.batch_size : (idx + 1) * self.batch_size])

        # Fetch only required tracks
        attrs = [self.label]
        q = MultiRangeQuery(self.momics, subrg)
        if self.data != "nucleotide":
            attrs.append(self.data)

        if self.silent:
            logging.disable(logging.WARNING)
        q.query_tracks(tracks=attrs)
        if self.silent:
            logging.disable(logging.WARNING)
        logging.disable(logging.NOTSET)

        # If input is a track, reshape and filter out NaN values
        if self.data in q.coverage.keys():  # type: ignore
            X = np.array(list(q.coverage[self.data].values()))  # type: ignore
            # filter = ~np.isnan(X).any(axis=1)
            # X = X[filter]
            X = np.nan_to_num(X, nan=0)
            sh = X.shape
            X = X.reshape(-1, sh[1], 1)

        # If input is the sequences, one-hot-encode the sequences and resize
        elif self.data == "nucleotide":
            q.query_sequence()
            seqs = list(q.seq["nucleotide"].values())  # type: ignore

            # One-hot-encode the sequences lists in seqs
            def one_hot_encode(seq) -> np.ndarray:
                seq = seq.upper()
                encoding_map = {"A": [1, 0, 0, 0], "T": [0, 1, 0, 0], "C": [0, 0, 1, 0], "G": [0, 0, 0, 1]}
                oha = np.zeros((len(seq), 4), dtype=int)
                for i, nucleotide in enumerate(seq):
                    oha[i] = encoding_map[nucleotide]

                return oha

            X = np.array([one_hot_encode(seq) for seq in seqs])
            sh = X.shape
            X = X.reshape(-1, sh[1], 4)
        else:
            raise ValueError("data must be a track label or 'nucleotide'")

        # Extract label and filter out NaN values
        out = np.array(list(q.coverage[self.label].values()))  # type: ignore
        # out = out[filter]
        out = np.nan_to_num(out, nan=0)

        # Recenter label if needed
        if self.label_size is not None:
            midpos = out.shape[1] // 2
            out = out[:, int(midpos - self.label_size / 2) : int(midpos + self.label_size / 2)]
            dim = self.label_size
        else:
            sh = out.shape
            dim = sh[1]

        Y = out.reshape(-1, dim, 1)

        return X, Y

    def __str__(self):
        return f"RangeDataLoader(start={self.start}, stop={self.stop}, batch_size={self.batch_size})"

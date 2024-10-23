from typing import Callable, Optional, Generator, Tuple

import numpy as np
import pyranges as pr
import logging

from .momics import Momics
from .logging import logger
from .momicsquery import MomicsQuery


class MomicsStreamer:
    """
    This class is implemented to efficiently query a `momics` repository by batches
    and extract any coverage data from it. The data streamer will iterate over ranges in batches
    and iteratively query a `momics`.

    For a tensorflow DataSet constructor, see `momics.dataset.MomicsDataset`.

    See Also:
        :class:`momics.dataset.MomicsDataset`

    Attributes:
        momics (Momics): a local `.momics` repository.
        ranges (dict): pr.PyRanges object.
        batch_size (int): the batch size
        features (list): list of track labels to query
        silent (bool): whether to suppress info messages
    """

    def __init__(
        self,
        momics: Momics,
        ranges: pr.PyRanges,
        batch_size: Optional[int] = None,
        features: Optional[int] = None,
        preprocess_func: Optional[Callable] = None,
        silent: bool = True,
    ) -> None:
        """Initialize the MomicsStreamer object.

        Args:
            momics (Momics): a Momics object
            ranges (dict): pr.PyRanges object.
            batch_size (int): the batch size
            features (list): list of track labels to query
            preprocess_func (Callable): a function to preprocess the queried data
            silent (bool): whether to suppress info messages
        """

        self.momics = momics
        self.ranges = ranges
        if batch_size is None:
            batch_size = len(ranges)
        self.batch_size = batch_size
        self.num_batches = (len(ranges) + batch_size - 1) // batch_size
        if features is not None:
            if not isinstance(features, list):
                features = [features]
            i = len(features)
            if "nucleotide" in features:
                i -= 1
                _ = momics.seq()  # Check that the momics object has a sequence
            if i > 0:  # Other features besides "nucleotide"
                tr = momics.tracks()  # Check that the momics object has the tracks
                for f in features:
                    if f == "nucleotide":
                        continue
                    if f not in list(tr["label"]):
                        raise ValueError(f"Features {f} not found in momics repository.")

        self.features = features
        self.silent = silent
        self.preprocess_func = preprocess_func if preprocess_func else self._default_preprocess
        self.batch_index = 0

    def query(self, batch_ranges) -> Tuple:
        """
        Query function to fetch data from a `momics` repo based on batch_ranges.

        Args:
            batch_ranges (pr.PyRanges): PyRanges object for a batch

        Returns:
            Queried coverage/sequence data
        """

        attrs = self.features
        i = len(attrs)
        res = {attr: None for attr in attrs}
        q = MomicsQuery(self.momics, batch_ranges)

        if self.silent:
            logging.disable(logging.WARNING)

        # Fetch seq if needed
        if "nucleotide" in attrs:
            i -= 1
            q.query_sequence()
            seqs = list(q.seq["nucleotide"].values())

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
            res["nucleotide"] = X.reshape(-1, sh[1], 4)

        # Fetch coverage tracks if needed
        if i > 0:
            attrs2 = [attr for attr in attrs if attr != "nucleotide"]
            q.query_tracks(tracks=attrs2)
            for attr in attrs2:
                out = np.array(list(q.coverage[attr].values()))
                sh = out.shape
                res[attr] = out.reshape(-1, sh[1], 1)

        if self.silent:
            logging.disable(logging.NOTSET)

        return tuple(res.values())

    def _default_preprocess(self, data):
        """
        Default preprocessing function that normalizes data.
        """
        return (data - np.mean(data, axis=0)) / np.std(data, axis=0)

    def generator(self) -> Generator:
        """
        Generator to yield batches of ranges and queried/preprocessed data.

        Yields:
            Tuple[pr.PyRanges, np.ndarray]: batch_ranges and preprocessed_data
        """
        self.batch_index = 0
        for i in range(0, len(self.ranges), self.batch_size):
            batch_ranges = pr.PyRanges(self.ranges.df.iloc[i : i + self.batch_size])
            queried_data = self.query(batch_ranges)
            # preprocessed_data = self.preprocess(queried_data)
            self.batch_index += 1
            yield queried_data

    def __iter__(self):
        return self.generator()

    def __next__(self):
        """Return the next batch or raise StopIteration."""
        if self.batch_index < self.num_batches:
            start = self.batch_index * self.batch_size
            end = min((self.batch_index + 1) * self.batch_size, len(self.ranges))
            batch_ranges = pr.PyRanges(self.ranges.df.iloc[start:end])
            queried_data = self.query(batch_ranges)
            self.batch_index += 1
            return queried_data
        else:
            raise StopIteration

    def __len__(self):
        return self.num_batches

    def reset(self):
        """Reset the iterator to allow re-iteration."""
        self.batch_index = 0

    def batch(self, batch_size: int):
        """
        Change the batch size for streaming data.

        Args:
            batch_size (int): The new size for batches.
        """
        if batch_size <= 0:
            raise ValueError("Batch size must be greater than zero.")

        if batch_size > len(self.ranges):
            batch_size = len(self.ranges)

        self.batch_size = batch_size
        self.num_batches = (len(self.ranges) + self.batch_size - 1) // self.batch_size
        self.reset()
        logger.info(f"Batch size updated to {self.batch_size}. Number of batches is now {self.num_batches}.")

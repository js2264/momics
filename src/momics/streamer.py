from typing import Callable, Optional, Generator, Union

import numpy as np
import pyranges as pr

from . import utils as mutils
from .momics import Momics
from .logging import logger
from .query import MomicsQuery


class MomicsStreamer:
    """
    This class is implemented to efficiently query a `momics` repository for a set of features by batch,
    and extract data from it. The data streamer will iterate over ranges in batches
    and iteratively query a `momics`. The streamer can also be used to preprocess the data before returning it.

    When iterating over the streamer by batch, each iteration will return a dictionary with the queried data.
    The keys of the dictionary are the queried features and the values are numpy arrays
    containing the data over the corresponding ranges.

    For a tensorflow DataSet constructor, see `momics.dataset.MomicsDataset`.

    See Also:
        :class:`momics.dataset.MomicsDataset`

    Attributes:
        momics (Momics): a local `.momics` repository.
        ranges (dict): pr.PyRanges object.
        batch_size (int): the batch size
        features (list): list of track labels to query
        silent (bool): whether to suppress info messages

    Example:
        >>> from momics import streamer as mms
        >>> from momics import momics as mmm
        >>> from pyranges import PyRanges
        >>>
        >>> # Create a MomicsStreamer object
        >>> repo = mmm.Momics("yeast.momics")
        >>> ranges = repo.bins(1000)
        >>> stream = mms.MomicsStreamer(repo, ranges, batch_size=1000, features=["nucleotide", "atac"])
        >>>
        >>> # Iterate over the streamer
        >>> for batch in stream:
        >>>     # Process the batch
        >>>     print(batch)
    """

    def __init__(
        self,
        momics: Momics,
        ranges: pr.PyRanges,
        batch_size: Optional[int] = None,
        features: Optional[Union[list, str]] = None,
        preprocess_func: Optional[Callable] = None,
        silent: bool = True,
    ) -> None:
        """Initialize the MomicsStreamer object.

        Args:
            momics (Momics): a Momics object
            ranges (dict): pr.PyRanges object
            batch_size (int): the batch size
            features (list): list of track labels to query
            preprocess_func (Callable): a function to preprocess the queried data
            silent (bool): whether to suppress info messages
        """

        self.momics = momics
        self.ranges = ranges

        # Set the batch size
        if batch_size is None:
            batch_size = len(ranges)

        self.batch_size = batch_size
        self.num_batches = (len(ranges) + batch_size - 1) // batch_size

        # Check features
        if features is not None:
            if isinstance(features, str):
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
        else:
            features = list(momics.tracks()["label"])
        self.features = features
        self.silent = silent
        self.preprocess_func = preprocess_func if preprocess_func else self._no_preprocess
        self.batch_index = 0

    def query(self, batch_ranges) -> dict:
        """
        Query function to fetch data from a `momics` repo based on batch_ranges.

        Args:
            batch_ranges (pr.PyRanges): PyRanges object for a batch

        Returns:
            Queried coverage/sequence data
        """

        attrs = self.features
        i = len(attrs)
        res: dict = {attr: None for attr in attrs}
        q = MomicsQuery(self.momics, batch_ranges)

        # Fetch seq if needed
        if "nucleotide" in attrs:
            i -= 1
            q.query_sequence()
            if q.seq is not None:
                seqs = list(q.seq["nucleotide"].values())
            else:
                raise ValueError("No sequence data found in the momics repository.")

            X = np.array([mutils.one_hot_encode(seq) for seq in seqs])
            sh = X.shape
            res["nucleotide"] = X.reshape(-1, sh[1], 4)

        # Fetch coverage tracks if needed
        if i > 0:
            attrs2 = [attr for attr in attrs if attr != "nucleotide"]
            q.query_tracks(tracks=attrs2)
            for attr in attrs2:
                if q.coverage is not None:
                    out = np.array(list(q.coverage[attr].values()))
                else:
                    raise ValueError(f"{attr} track data found in the momics repository.")

                sh = out.shape
                res[attr] = out.reshape(-1, sh[1], 1)

        return res

    def _zscore(self, data):
        """
        Z-score preprocessing function to normalize data
        """
        return (data - np.mean(data, axis=0)) / np.std(data, axis=0)

    def _no_preprocess(self, data):
        """
        No preprocessing function
        """
        return data

    def __iter__(self):
        self.batch_index = 0
        return self

    def __next__(self):
        """Return the next batch or raise StopIteration."""
        if self.batch_index >= self.num_batches:
            raise StopIteration
        start = self.batch_index * self.batch_size
        end = min((self.batch_index + 1) * self.batch_size, len(self.ranges))
        batch_ranges = pr.PyRanges(self.ranges.df.iloc[start:end])
        self.batch_index += 1
        return self.query(batch_ranges)

    def __len__(self):
        return self.num_batches

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
        self.batch_index = 0
        logger.debug(f"Batch size updated to {self.batch_size}. Number of batches is now {self.num_batches}.")

    def generator(self) -> Generator:
        """
        Generator to yield batches of data.
        """
        for i in range(0, len(self.ranges), self.batch_size):
            batch_ranges = pr.PyRanges(self.ranges.df.iloc[i : i + self.batch_size])
            yield self.query(batch_ranges)

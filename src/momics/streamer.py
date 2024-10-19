from typing import Callable, Optional

import numpy as np
import pyranges as pr
import logging

from .momics import Momics
from .multirangequery import MultiRangeQuery


class MomicsStreamer:
    """
    This class is implemented to efficiently query a `momics` repository by batches
    and extract any coverage data from it. The data streamer will iterate over ranges in batches
    and iteratively query a `momics`.

    For a tensorflow DataSet constructor, see `momics.dataset.MomicsDataset`.

    See Also
    --------
        `momics.dataset.MomicsDataset`

    Attributes
    ----------
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
        silent: bool = False,
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
            i = 0
            if not isinstance(features, list):
                features = [features]
            if "nucleotide" in features:
                features.remove("nucleotide")
                i += 1
                _ = momics.seq()
            if len(features) > i:
                tr = momics.tracks()
                for f in features:
                    if f not in list(tr["label"]):
                        raise ValueError(f"Features {f} not found in momics repository.")

        if i > 0:
            features.insert(0, "nucleotide")
        self.features = features
        self.silent = silent
        self.preprocess_func = preprocess_func if preprocess_func else self._default_preprocess

    def query(self, batch_ranges):
        """
        Query function to fetch data from a `momics` repo based on batch_ranges.

        Args:
            batch_ranges (pr.PyRanges): PyRanges object for a batch

        Returns:
            Queried coverage/sequence data
        """

        attrs = self.features
        res = {attr: None for attr in attrs}
        q = MultiRangeQuery(self.momics, batch_ranges)

        if self.silent:
            logging.disable(logging.WARNING)

        # Fetch seq if needed
        if "nucleotide" in attrs:
            attrs.remove("nucleotide")
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
        if len(attrs) > 0:
            q.query_tracks(tracks=attrs)
            for attr in attrs:
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

    def generator(self):
        """
        Generator to yield batches of ranges and queried/preprocessed data.

        Yields:
            Tuple[pr.PyRanges, np.ndarray]: batch_ranges and preprocessed_data
        """
        for i in range(0, len(self.ranges), self.batch_size):
            batch_ranges = pr.PyRanges(self.ranges.df.iloc[i : i + self.batch_size])
            queried_data = self.query(batch_ranges)
            # preprocessed_data = self.preprocess(queried_data)
            yield queried_data

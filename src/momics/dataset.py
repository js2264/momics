from typing import Callable, Generator, Optional, Tuple

import numpy as np
import pyranges as pr
import logging
import tensorflow as tf

from .momics import Momics
from .streamer import MomicsStreamer
from .multirangequery import MultiRangeQuery


class MomicsDataset(tf.data.Dataset):
    """
    This class is implemented to train deep learning models, where the
    input data (features) is a track or a sequence and the labeled data (target)
    is another track. The data loader will iterate over the ranges in batches
    and extract the features and target for each range. It is a subclass of
    `tf.data.DataSet` and can be used as a generator for a `tf.keras.Model`.

    For a more basic generator to stream a `momics` by batches of ranges,
    see `momics.streamer.MomicsStreamer`.

    See Also
    --------
        `momics.streamer.MomicsStreamer`
    """

    def __new__(
        cls,
        momics: Momics,
        ranges: pr.PyRanges,
        features: str,
        target: str,
        target_size: Optional[int] = None,
        batch_size: Optional[int] = None,
        preprocess_func: Optional[Callable] = None,
        shuffle_buffer_size: int = 10000,
        prefetch_buffer_size: Optional[int] = tf.data.experimental.AUTOTUNE,
        silent: bool = False,
    ) -> tf.data.Dataset:
        """Create the MomicsDataset object.

        Args:
            momics (Momics): a Momics object
            ranges (pr.PyRanges): pr.PyRanges object
            features (str): the name of the track to use for input data
            target_size (int): To which width should the target be centered
            target (str): the name of the track to use for output data
            batch_size (int): the batch size
            preprocess_func (Callable): a function to preprocess the queried data
            shuffle_buffer_size (int): the size of the shuffle buffer. Pass 0 to disable shuffling
            prefetch_buffer_size (int): the size of the prefetch buffer
            silent (bool): whether to suppress info messages
        """

        # Check that all ranges have the same width
        df = ranges.df
        widths = df.End - df.Start
        if len(set(widths)) != 1:
            raise ValueError("All ranges must have the same width")
        w = int(widths[0])

        # Check that the target size is smaller than the features width
        if target_size is not None and target_size > w:
            raise ValueError("Target size must be smaller than the features width.")

        # Encapsulate MomicsStreamer logic
        def generator() -> Generator:
            streamer = MomicsStreamer(
                momics, ranges, batch_size, features=[features, target], preprocess_func=preprocess_func, silent=silent
            ).generator()
            for features_data, out in streamer:

                # Adjust the output if target_size is provided
                if target_size:
                    center = out.shape[1] // 2
                    label_data = out[:, int(center - target_size // 2) : int(center + target_size // 2)]
                else:
                    label_data = out

                yield features_data, label_data

        # Example output signature (modify based on your actual data shapes)
        feature_shape = (None, w, 4 if features == "nucleotide" else 1)
        label_shape = (None, target_size if target_size else w, 4 if features == "nucleotide" else 1)

        dataset = tf.data.Dataset.from_generator(
            generator,
            output_signature=(
                tf.TensorSpec(shape=feature_shape, dtype=tf.float32),
                tf.TensorSpec(shape=label_shape, dtype=tf.float32),
            ),
        )

        # Add shuffling and prefetching
        if shuffle_buffer_size > 0:
            shuffle_buffer_size = min(shuffle_buffer_size, batch_size)
            dataset = dataset.shuffle(shuffle_buffer_size)

        dataset = dataset.prefetch(buffer_size=prefetch_buffer_size)

        return dataset


class RangeDataLoader(tf.keras.utils.Sequence):
    """
    This class is implemented to train deep learning models, where the
    input data (features) is a track or a sequence and the labeled data (target)
    is another track. The data loader will iterate over the ranges in batches
    and extract the features and target for each range. It is a subclass of
    `tf.data.DataSet` and can be used as a generator for a `tf.keras.Model`.

    For a more basic generator to stream a `momics` by batches of ranges,
    see `momics.streamer.MomicsStreamer`.

    See Also
    --------
        `momics.streamer.MomicsStreamer`

    Attributes
    ----------
    momics (Momics): a local `.momics` repository.
    ranges (dict): pr.PyRanges object.
    features (str): the name of the track to use for input data
    target (str): the name of the track to use for output data
    target_size (int): To which width should the target be centered
    """

    def __init__(
        self,
        momics: Momics,
        ranges: pr.PyRanges,
        features: str,
        target: str,
        target_size: Optional[int] = None,
        batch_size: Optional[int] = None,
        silent: bool = False,
    ) -> None:
        """Initialize the RangeDataLoader object.

        Args:
            momics (Momics): a Momics object
            ranges (pr.PyRanges): pr.PyRanges object
            features (str): the name of the track to use for input data
            target_size (int): To which width should the target be centered
            target (str): the name of the track to use for output data
            batch_size (int): the batch size
            silent (bool): whether to suppress info messages
        """

        # Check that all ranges have the same width
        df = ranges.df
        widths = df.End - df.Start
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
        if features == "nucleotide":
            _ = momics.seq()

        if features not in list(tr["label"]) and features != "nucleotide":
            raise ValueError(f"Features {features} not found in momics repository.")
        if target not in list(tr["label"]):
            raise ValueError(f"Target {target} not found in momics repository.")

        self.features = features
        self.target = target

        if target_size is not None and target_size > int(widths[0]):
            raise ValueError("Target size must be smaller than the features width.")
        self.target_size = target_size

    def __len__(self) -> int:
        return int(np.ceil(len(self.ranges) / self.batch_size))

    def __getitem__(self, idx) -> Tuple[np.ndarray, np.ndarray]:
        subrg = pr.PyRanges(self.ranges.df[idx * self.batch_size : (idx + 1) * self.batch_size])

        # Fetch only required tracks
        attrs = [self.target]
        q = MultiRangeQuery(self.momics, subrg)
        if self.features != "nucleotide":
            attrs.append(self.features)

        if self.silent:
            logging.disable(logging.WARNING)
        q.query_tracks(tracks=attrs)
        if self.silent:
            logging.disable(logging.WARNING)
        logging.disable(logging.NOTSET)

        # If input is a track, reshape and filter out NaN values
        if self.features in q.coverage.keys():  # type: ignore
            X = np.array(list(q.coverage[self.features].values()))  # type: ignore
            # filter = ~np.isnan(X).any(axis=1)
            # X = X[filter]
            X = np.nan_to_num(X, nan=0)
            sh = X.shape
            X = X.reshape(-1, sh[1], 1)

        # If input is the sequences, one-hot-encode the sequences and resize
        elif self.features == "nucleotide":
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
            raise ValueError("features must be a track label or 'nucleotide'")

        # Extract label and filter out NaN values
        out = np.array(list(q.coverage[self.target].values()))  # type: ignore
        # out = out[filter]
        out = np.nan_to_num(out, nan=0)

        # Recenter label if needed
        if self.target_size is not None:
            midpos = out.shape[1] // 2
            out = out[:, int(midpos - self.target_size / 2) : int(midpos + self.target_size / 2)]
            dim = self.target_size
        else:
            sh = out.shape
            dim = sh[1]

        Y = out.reshape(-1, dim, 1)

        return X, Y

    def __str__(self):
        return f"RangeDataLoader(start={self.start}, stop={self.stop}, batch_size={self.batch_size})"

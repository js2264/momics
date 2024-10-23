from typing import Callable, Generator, Optional

import pyranges as pr
import tensorflow as tf

from .momics import Momics
from .streamer import MomicsStreamer


class MomicsDataset(tf.data.Dataset):
    """
    This class is implemented to train deep learning models, where the
    input data (features) is a track or a sequence and the labeled data (target)
    is another track. The data loader will iterate over the ranges in batches
    and extract the features and target for each range. It is a subclass of
    `tf.data.DataSet` and can be used as a generator for a `tf.keras.Model`.

    For a more basic generator to stream a `momics` by batches of ranges,
    see `momics.streamer.MomicsStreamer`.

    See Also:
        :class:`momics.streamer.MomicsStreamer`
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
        silent: bool = True,
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
            )
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
        label_shape = (None, target_size if target_size else w, 4 if target == "nucleotide" else 1)

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

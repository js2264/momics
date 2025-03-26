from typing import Callable, Optional

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
        repo: Momics,
        ranges: pr.PyRanges,
        features: str,
        target: str,
        target_size: Optional[int] = None,
        batch_size: Optional[int] = None,
        preprocess_func: Optional[Callable] = None,
        silent: bool = True,
        cache: bool = False,
    ) -> tf.data.Dataset:
        """Create the MomicsDataset object.

        Args:
            repo (Momics): a Momics object
            ranges (pr.PyRanges): pr.PyRanges object
            features (str): the name of the track to use for input data
            target (str): the name of the track to use for output data
            target_size (int): To which width should the target be centered
            batch_size (int): the batch size
            preprocess_func (Callable): a function to preprocess the queried data
            silent (bool): whether to suppress info messages
            cache (bool): whether to cache the dataset
        """

        # Check that all ranges have the same width
        widths = ranges.End - ranges.Start
        if len(set(widths)) != 1:
            raise ValueError("All ranges must have the same width")
        features_size = int(widths.iloc[0])

        # Check that the target size is smaller than the features width
        if target_size is None:
            target_size = features_size
        if target_size is not None and target_size > features_size:
            raise ValueError(f"Target size must be smaller than or equal to the feature size {features_size}.")

        ranges_target = ranges.copy()
        if target_size != features_size:
            ranges_target.Start = ranges_target.Start + features_size // 2 - target_size // 2
            ranges_target.End = ranges_target.Start + target_size

        # Define generator for features data
        x_streamer = MomicsStreamer(repo, ranges, batch_size, features=[features], preprocess_func=preprocess_func, silent=silent)
        x_gen = x_streamer.generator
        if features == "nucleotide":
            out = tf.TensorSpec(shape=(None, features_size, 4), dtype=tf.int32)
        else:
            out = tf.TensorSpec(shape=(None, features_size, 1), dtype=tf.float32)

        x_dataset = tf.data.Dataset.from_generator(x_gen, output_signature=(out,))

        # Define generator for target data
        y_streamer = MomicsStreamer(
            repo, ranges_target, batch_size, features=[target], preprocess_func=preprocess_func, silent=silent
        )
        y_gen = y_streamer.generator
        if target == "nucleotide":
            out = tf.TensorSpec(shape=(None, target_size, 4), dtype=tf.int32)
        else:
            out = tf.TensorSpec(shape=(None, target_size, 1), dtype=tf.float32)

        y_dataset = tf.data.Dataset.from_generator(y_gen, output_signature=(out,))

        # Combine features and target datasets
        xy_ds = tf.data.Dataset.zip((x_dataset, y_dataset))

        # Add caching option
        if cache:
            xy_ds = xy_ds.cache()

        return xy_ds

    def __len__(self):
        return super().__len__()

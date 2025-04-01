from typing import Callable, Optional, Union

import pyranges as pr
import tensorflow as tf

from .momics import Momics
from .streamer import MomicsStreamer


class MomicsDataset(tf.data.Dataset):
    def __new__(
        cls,
        repo: Momics,
        ranges: pr.PyRanges,
        features: Union[list, str],
        target: Optional[Union[list, str]] = None,
        target_size: Optional[int] = None,
        batch_size: Optional[int] = None,
        preprocess_func: Optional[Callable] = None,
        silent: bool = True,
        cache: bool = False,
    ) -> tf.data.Dataset:
        """
        This class is implemented to train deep learning models, where the
        input data (features) is a track or a sequence and the labeled data (target)
        is another track. The data loader will iterate over the ranges in batches
        and extract the features and target for each range. It is a subclass of
        `tf.data.DataSet` and can be used as a generator for a `tf.keras.Model`.

        For a more basic generator to stream a `momics` by batches of ranges,
        see `momics.streamer.MomicsStreamer`.

        Args:
            repo (Momics): a Momics object
            ranges (pr.PyRanges): pr.PyRanges object
            features (str): the name of the track to use for input data (can also be "nucleotide")
            target (str): the name of the track to use for output data
            target_size (int): To which width should the target be centered
            batch_size (int): the batch size
            preprocess_func (Callable): a function to preprocess the queried data
            silent (bool): whether to suppress info messages
            cache (bool): whether to cache the dataset

        Returns:
            tf.data.Dataset: a TensorFlow dataset object. The object is a generator
            that yields batches of data.

            The yielded data is a nested tuple ((features1, features2, ...), (target1, target2, ...))
            where features1, features2, ... are the features of the input data and
            target1, target2, ... are the targets of the output data:

            - The shape of the yielded input data is (batch_size, features_size, 4) for nucleotide data and
            (batch_size, features_size, 1) for other data.
            - The shape of the yielded output data is (batch_size, target_size, 4) for nucleotide data and
            (batch_size, target_size, 1) for other data.
            - If target is None, the yielded data is a tuple of (features1, features2, ...).

        See Also:
            :class:`momics.streamer.MomicsStreamer`
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
        xds = []
        if isinstance(features, str):
            features = [features]
        for ft in features:
            x_streamer = MomicsStreamer(repo, ranges, batch_size, features=[ft], preprocess_func=preprocess_func, silent=silent)
            x_gen = x_streamer.generator
            if ft == "nucleotide":
                out = tf.TensorSpec(shape=(None, features_size, 4), dtype=tf.int32)
            else:
                out = tf.TensorSpec(shape=(None, features_size, 1), dtype=tf.float32)
            xds.append(tf.data.Dataset.from_generator(x_gen, output_signature=(out,)))
        if len(xds) == 1:
            x_dataset = xds[0]
        else:
            x_dataset = tf.data.Dataset.zip(tuple(xds))

        # Define generator for target data
        if target is not None:
            if isinstance(target, str):
                target = [target]
            yds = []
            for tg in target:
                y_streamer = MomicsStreamer(
                    repo, ranges_target, batch_size, features=[tg], preprocess_func=preprocess_func, silent=silent
                )
                y_gen = y_streamer.generator
                if tg == "nucleotide":
                    out = tf.TensorSpec(shape=(None, target_size, 4), dtype=tf.int32)
                else:
                    out = tf.TensorSpec(shape=(None, target_size, 1), dtype=tf.float32)

                yds.append(tf.data.Dataset.from_generator(y_gen, output_signature=(out,)))

            if len(yds) == 1:
                y_dataset = yds[0]
            else:
                y_dataset = tf.data.Dataset.zip(tuple(yds))

            xy_ds = tf.data.Dataset.zip((x_dataset, y_dataset))

        else:
            xy_ds = x_dataset

        # Add caching option
        if cache:
            xy_ds = xy_ds.cache()

        # Calculate the number of batches in the dataset
        if batch_size is not None:
            n_samples = len(ranges)
            n_batches = (n_samples + batch_size - 1) // batch_size
            xy_ds = xy_ds.apply(tf.data.experimental.assert_cardinality(n_batches))

        return xy_ds

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
        Create a TensorFlow dataset for training deep learning models on genomic data.

        This dataset generator is designed for training deep learning models on genomic data,
        supporting multiple input features (tracks or nucleotide sequences) and multiple
        output targets. It can be used to train models with any number of labels and heads.

        The dataset automatically handles batching, centering of target regions, and
        preprocessing of input data. It seamlessly integrates with TensorFlow's training
        pipelines as it is a subclass of `tf.data.Dataset`.

        Parameters
        ----------
        repo : Momics
            A Momics repository object containing the genomic data
        ranges : pr.PyRanges
            PyRanges object specifying the genomic regions to extract data from.
            All ranges must have the same width.
        features : Union[list, str]
            Track name(s) to use as input features. Can be a single string or a list of strings.
            Use "nucleotide" to include genomic sequence data.
        target : Optional[Union[list, str]], default=None
            Track name(s) to use as prediction targets. Can be a single string or a list of strings.
            If None, the dataset will only yield feature data.
        target_size : Optional[int], default=None
            Width to which the target should be centered. If None, uses the same size as features.
            Must be less than or equal to the feature size.
        batch_size : Optional[int], default=None
            Number of samples per batch. If None, the entire dataset is returned as one batch.
        preprocess_func : Optional[Callable], default=None
            Function to preprocess the queried data before yielding.
            Should accept the raw data and return the processed data.
        silent : bool, default=True
            Whether to suppress information messages during dataset creation.
        cache : bool, default=False
            Whether to cache the dataset in memory after the first iteration.

        Returns
        -------
        tf.data.Dataset
            A TensorFlow dataset that yields batches of data in the following format:

            - With targets: ({feature_dict}, {target_dict})
              Where {feature_dict} and {target_dict} are dictionaries mapping track names to data arrays
            - Without targets: {feature_dict}
              Where {feature_dict} is a dictionary mapping track names to data arrays

            Data shapes:
            - Nucleotide data: (batch_size, sequence_length, 4) with dtype=tf.int32
            - Track data: (batch_size, sequence_length, 1) with dtype=tf.float32

        Examples
        --------
        ```python
        # Create a dataset with one feature and one target
        dataset = MomicsDataset(
            repo=my_repo,
            ranges=my_ranges,
            features="ATAC-seq",
            target="ChIP-seq",
            batch_size=32
        )

        # Create a dataset with multiple features and targets
        dataset = MomicsDataset(
            repo=my_repo,
            ranges=my_ranges,
            features=["ATAC-seq", "nucleotide", "methylation"],
            target=["ChIP-seq", "expression"],
            batch_size=32,
            target_size=128,
            cache=True
        )

        # Use with Keras model
        model.fit(dataset, epochs=10)
        ```

        See Also
        --------
        momics.streamer.MomicsStreamer : For a more basic generator to stream data by batches
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

        # Adjust the target ranges to center of features
        ranges_target = ranges.copy()
        if target_size != features_size:
            ranges_target.Start = ranges_target.Start + features_size // 2 - target_size // 2
            ranges_target.End = ranges_target.Start + target_size

        # Check if features/targets are lists of strings
        if isinstance(features, str):
            features = [features]
        if isinstance(target, str):
            target = [target]

        # Define list of streamers for features data
        xstreamers = {}
        for ft in features:
            xstreamers[ft] = MomicsStreamer(
                repo, ranges, batch_size, features=[ft], preprocess_func=preprocess_func, silent=silent
            )

        # Define list of streamers for target data
        ystreamers = {}
        if target is not None:
            for tg in target:
                ystreamers[tg] = MomicsStreamer(
                    repo, ranges_target, batch_size, features=[tg], preprocess_func=preprocess_func, silent=silent
                )

        # Callable for combined generator
        def combined_generator():
            """Generate combined features and targets."""
            feature_generators = {ft: xstreamers[ft].generator() for ft in features}
            if target is not None:
                target_generators = {tg: ystreamers[tg].generator() for tg in target}
            else:
                target_generators = {}

            primary_feature = features[0]
            primary_gen = feature_generators[primary_feature]

            for x_batch in primary_gen:
                feature_batches = {}
                target_batches = {}

                feature_batches[primary_feature] = x_batch[primary_feature]  # Store primary feature data
                for ft in features:
                    if ft == primary_feature:
                        continue  # Skip primary feature (already processed)

                    feature_batch = next(feature_generators[ft])
                    feature_batches[ft] = feature_batch[ft]

                if target is not None:
                    for tg in target:
                        target_batch = next(target_generators[tg])
                        target_batches[tg] = target_batch[tg]

                yield (feature_batches, target_batches)

        # Define output signatures
        xsigs = {
            ft: tf.TensorSpec(
                shape=(None, features_size, 4) if ft == "nucleotide" else (None, features_size, 1),
                dtype=tf.int32 if ft == "nucleotide" else tf.float32,
                name=ft,
            )
            for ft in features
        }
        if target is not None:
            ysigs = {
                tg: tf.TensorSpec(
                    shape=(None, target_size, 4) if tg == "nucleotide" else (None, target_size, 1),
                    dtype=tf.int32 if tg == "nucleotide" else tf.float32,
                    name=tg,
                )
                for tg in target
            }
        else:
            ysigs = {}

        out_sigs = (xsigs, ysigs)

        # Create the final dataset from the generator
        xy_ds = tf.data.Dataset.from_generator(combined_generator, output_signature=out_sigs)

        # Replace xy_ds if target is None
        if target is None:
            xy_ds = xy_ds.map(lambda x, _: x)

        # Add caching option
        if cache:
            xy_ds = xy_ds.cache()

        # Calculate the number of batches in the dataset
        if batch_size is not None:
            n_samples = len(ranges)
            n_batches = (n_samples + batch_size - 1) // batch_size
            xy_ds = xy_ds.apply(tf.data.experimental.assert_cardinality(n_batches))

        return xy_ds

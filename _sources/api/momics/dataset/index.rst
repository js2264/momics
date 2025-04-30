momics.dataset
==============

.. py:module:: momics.dataset


Classes
-------

.. autoapisummary::

   momics.dataset.MomicsDataset


Module Contents
---------------

.. py:class:: MomicsDataset(variant_tensor)

   Bases: :py:obj:`tensorflow.data.Dataset`


   Represents a potentially large set of elements.

   The `tf.data.Dataset` API supports writing descriptive and efficient input
   pipelines. `Dataset` usage follows a common pattern:

   1. Create a source dataset from your input data.
   2. Apply dataset transformations to preprocess the data.
   3. Iterate over the dataset and process the elements.

   Iteration happens in a streaming fashion, so the full dataset does not need to
   fit into memory.

   Source Datasets:

   The simplest way to create a dataset is to create it from a python `list`:

   >>> dataset = tf.data.Dataset.from_tensor_slices([1, 2, 3])
   >>> for element in dataset:
   ...   print(element)
   tf.Tensor(1, shape=(), dtype=int32)
   tf.Tensor(2, shape=(), dtype=int32)
   tf.Tensor(3, shape=(), dtype=int32)

   To process lines from files, use `tf.data.TextLineDataset`:

   >>> dataset = tf.data.TextLineDataset(["file1.txt", "file2.txt"])

   To process records written in the `TFRecord` format, use `TFRecordDataset`:

   >>> dataset = tf.data.TFRecordDataset(["file1.tfrecords", "file2.tfrecords"])

   To create a dataset of all files matching a pattern, use
   `tf.data.Dataset.list_files`:

   ```python
   dataset = tf.data.Dataset.list_files("/path/*.txt")
   ```

   See `tf.data.FixedLengthRecordDataset` and `tf.data.Dataset.from_generator`
   for more ways to create datasets.

   Transformations:

   Once you have a dataset, you can apply transformations to prepare the data for
   your model:

   >>> dataset = tf.data.Dataset.from_tensor_slices([1, 2, 3])
   >>> dataset = dataset.map(lambda x: x*2)
   >>> [a.item() for a in dataset.as_numpy_iterator()]
   [2, 4, 6]

   Common Terms:

   **Element**: A single output from calling `next()` on a dataset iterator.
     Elements may be nested structures containing multiple components. For
     example, the element `(1, (3, "apple"))` has one tuple nested in another
     tuple. The components are `1`, `3`, and `"apple"`.

   **Component**: The leaf in the nested structure of an element.

   Supported types:

   Elements can be nested structures of tuples, named tuples, and dictionaries.
   Note that Python lists are *not* treated as nested structures of components.
   Instead, lists are converted to tensors and treated as components. For
   example, the element `(1, [1, 2, 3])` has only two components; the tensor `1`
   and the tensor `[1, 2, 3]`. Element components can be of any type
   representable by `tf.TypeSpec`, including `tf.Tensor`, `tf.data.Dataset`,
   `tf.sparse.SparseTensor`, `tf.RaggedTensor`, and `tf.TensorArray`.

   ```python
   a = 1 # Integer element
   b = 2.0 # Float element
   c = (1, 2) # Tuple element with 2 components
   d = {"a": (2, 2), "b": 3} # Dict element with 3 components
   Point = collections.namedtuple("Point", ["x", "y"])
   e = Point(1, 2) # Named tuple
   f = tf.data.Dataset.range(10) # Dataset element
   ```

   For more information,
   read [this guide](https://www.tensorflow.org/guide/data).



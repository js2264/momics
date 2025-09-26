momics.attribution
==================

.. py:module:: momics.attribution


Functions
---------

.. autoapisummary::

   momics.attribution.analyze_mutation_impacts
   momics.attribution.attribution
   momics.attribution.get_ISM
   momics.attribution.get_saliency
   momics.attribution.mutate_sequence


Module Contents
---------------

.. py:function:: analyze_mutation_impacts(original_pred, mutated_preds, mutation_info, viewpoint_width = 200)

   Analyze the impact of each mutation on ATAC predictions.

   :param original_pred: Score predicted from the original sequence
   :param mutated_preds: Predictions for mutated sequences
   :param mutation_info: List of tuples (position, original_nuc, mutated_nuc)
   :param viewpoint_width: Size of window around mutation site to average predictions

   :returns: Dictionary with position as key and nucleotide impacts as values
   :rtype: ism_pos


.. py:function:: attribution(repo, model, track_name, chromosome, centerpoint, viewpoint_width = 0, batch = 0, skip_ism = False)

   Perform in-silico saturation mutagenesis (ISM) and calculate saliency
   of a given head over a given genomic region.

   :param repo: Momics repository object
   :param model: Trained TensorFlow/Keras model for prediction
   :param track_name: Name of the track to analyze (e.g., "bw2")
   :param chromosome: Chromosome name (e.g., "chr1")
   :param centerpoint: Center point of the region of interest (int)
   :param viewpoint_width: Width of the region around the centerpoint to analyze (int)
   :param batch: Batch size for model prediction. If 0, no batching is applied (int)
   :param skip_ism: If True, skip ISM computation and return None for ISM-related outputs (bool)

   :returns: One-hot encoded sequence of the region of interest
             data: Experimental data for the region of interest
             pred: Model predictions for the region of interest
             saliency: gradient x input of the model
             mutation_info: List of tuples (position, original_nuc, mutated_nuc)
             ism: a wx4 array, where w is the number of positions mutated, and the 4 columns are effects of mutating to A, T, G, C
             ism_score: the aggregated effect score for each position (average of the 3 mutations that are not the original nucleotide)
             ism_data: ism_score multiplied by the original experimental scores
   :rtype: seq


.. py:function:: get_ISM(model, seq, viewpoint_width=0, batch=0)

   Perform in-silico saturation mutagenesis (ISM) on the input sequence using the provided model.

   :param model: Trained TensorFlow/Keras model for prediction
   :param seq: One-hot encoded sequence, shape (1, x, 4)
   :param viewpoint_width: Width of the region around the centerpoint to analyze. If 0, uses model output width.
   :param batch: Batch size for model prediction. If 0, no batching is applied

   :returns: Dictionary with position as key and nucleotide impacts as values
             ism: Array of shape (w, 3) with nucleotide effects for each position
             ism_score: List of average ISM score for each position
   :rtype: ism_pos


.. py:function:: get_saliency(model, seq)

.. py:function:: mutate_sequence(seq, start, end)

   Generate all possible single nucleotide mutations in the sequence from start to end.
   :param seq: One-hot encoded sequence, shape (x, 4)
   :param start: Start position for mutations (inclusive)
   :param end: End position for mutations (exclusive)

   :returns: Array of mutated sequences, shape (800, x, 4)
             mutation_info: List of tuples (position, original_nuc, mutated_nuc)
   :rtype: mutated_sequences



import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import momics.query as mmq
import momics.utils as mutils


DEFAULT_ISM_CMAP = plt.get_cmap("afmhot_r")


def mutate_sequence(seq: np.ndarray, start: int, end: int) -> tuple[np.ndarray, list]:
    """
    Generate all possible single nucleotide mutations in the sequence from start to end.
    Args:
        seq: One-hot encoded sequence, shape (x, 5)
        start: Start position for mutations (inclusive)
        end: End position for mutations (exclusive)
    Returns:
        mutated_sequences: Array of mutated sequences, shape (800, x, 5)
        mutation_info: List of tuples (position, original_nuc, mutated_nuc)
    """
    mutated_sequences = []
    mutation_info = []
    nuc_names = ["N", "A", "T", "G", "C"]

    for pos in range(start, end):
        original_nuc_idx = np.argmax(seq[pos, :])
        for new_nuc_idx in range(1, 5):
            mutated_seq = np.copy(seq).reshape(1, seq.shape[0], seq.shape[1])
            mutated_seq[0, pos, :] = 0
            mutated_seq[0, pos, new_nuc_idx] = 1
            mutated_sequences.append(mutated_seq[0])
            mutation_info.append((pos, nuc_names[original_nuc_idx], nuc_names[new_nuc_idx]))

    return np.array(mutated_sequences), mutation_info


def analyze_mutation_impacts(
    original_pred: np.ndarray, mutated_preds: np.ndarray, mutation_info: list, viewpoint_width: int = 200
) -> dict:
    """
    Analyze the impact of each mutation on ATAC predictions.

    Args:
        original_pred: Score predicted from the original sequence
        mutated_preds: Predictions for mutated sequences
        mutation_info: List of tuples (position, original_nuc, mutated_nuc)
        viewpoint_width: Size of window around mutation site to average predictions

    Returns:
        ism_pos: Dictionary with position as key and nucleotide impacts as values
    """
    pred_center = original_pred.shape[0] // 2
    s = pred_center - viewpoint_width // 2
    e = pred_center + viewpoint_width // 2
    ism_pos = {}

    # Group mutations by position
    for i, (pos, orig_nuc, mut_nuc) in enumerate(mutation_info):
        if pos not in ism_pos:
            ism_pos[pos] = {"original_nuc": orig_nuc, "A": 0, "T": 0, "G": 0, "C": 0, "effect": 0}

        ism_pos[pos][mut_nuc] = np.abs(np.mean(mutated_preds[i][s:e]) - np.mean(original_pred[s:e]))

    for pos, x in ism_pos.items():
        orig_nuc = x["original_nuc"]
        effects = [x[nuc] for nuc in ["A", "T", "G", "C"] if nuc != orig_nuc]
        ism_pos[pos]["effect"] = np.mean(effects)

    return ism_pos


def plot_ISM_heatmap(
    ism_pos,
    figsize=(50, 2),
    cmap=DEFAULT_ISM_CMAP,
    figname="tmp.pdf",
):
    """
    Plot a heatmap of mutation impacts.
    Args:
        ism_pos: Dictionary with position as key and nucleotide impacts as values
        figsize: Size of the figure
        cmap: Colormap for the heatmap
        figname: Filename to save the figure
    """
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=figsize)
    positions = sorted(ism_pos.keys())
    positions_labs = [next(iter(ism_pos[pos].values())) for pos in positions]
    impacts = np.array([list(ism_pos[pos].values())[1:] for pos in positions])
    impacts0 = np.clip(impacts, np.quantile(impacts, 0.01), np.quantile(impacts, 0.99))

    cax = ax1.imshow(impacts0.T, aspect="auto", cmap=cmap, interpolation="nearest")
    ax1.set_xticks(np.arange(len(positions)))
    ax1.set_yticks(np.arange(4))
    ax1.set_yticklabels(["A", "T", "G", "C"])
    ax1.set_ylabel("Nucleotide")
    ax1.set_title("Impact of Mutations on ATAC Predictions")
    fig.colorbar(cax, ax=ax1, orientation="vertical", label="Impact Score")

    # Plot a second heatmap, below the first one
    impact_merged = np.sum(impacts, axis=1)
    cax = ax2.imshow(impact_merged.reshape(1, -1), aspect="auto", cmap=plt.get_cmap("afmhot_r"), interpolation="nearest")
    ax2.set_xticks(np.arange(len(positions)))
    ax2.set_xticklabels(positions_labs)
    ax2.set_yticks(np.arange(1))
    ax2.set_yticklabels(["N"])
    ax2.set_xlabel("Position")
    ax2.set_ylabel("Nucleotide")
    fig.colorbar(cax, ax=ax2, orientation="vertical", label="Impact Score")

    plt.savefig(figname, dpi=300, bbox_inches="tight")


def get_ISM(model, seq, viewpoint_width=0, batch=0):
    """
    Perform in-silico saturation mutagenesis (ISM) on the input sequence using the provided model.

    Args:
        model: Trained TensorFlow/Keras model for prediction
        seq: One-hot encoded sequence, shape (1, x, 5)
        viewpoint_width: Width of the region around the centerpoint to analyze. If 0, uses model output width.
        batch: Batch size for model prediction. If 0, no batching is applied

    Returns:
        ism_pos: Dictionary with position as key and nucleotide impacts as values
        ism: Array of shape (w, 4) with nucleotide effects for each position
        ism_score: List of average ISM score for each position
    """

    ## Find coordinates
    in_width = model.input_shape[1]
    out_width = model.output_shape[1]
    if viewpoint_width == 0:
        viewpoint_width = out_width

    # Mutate sequence
    seq_mutated, mutation_info = mutate_sequence(
        seq[0], in_width // 2 - viewpoint_width // 2, in_width // 2 + viewpoint_width // 2
    )

    # Run predictions on mutated sequences
    if batch > 0:
        seq_mutated_ds = tf.data.Dataset.from_generator(
            lambda: iter(seq_mutated), output_signature=tf.TensorSpec(shape=(in_width, 5), dtype=tf.float32)
        ).batch(batch)
    else:
        seq_mutated_ds = seq_mutated

    mutated_preds = model.predict(seq_mutated_ds, verbose=0)

    # Compute per-nucleotide attribution
    original_pred = model.predict(seq, verbose=0)
    ism_pos = analyze_mutation_impacts(original_pred[0], mutated_preds, mutation_info, viewpoint_width)

    # Format all outputs
    l_ism = []
    for pos in sorted(ism_pos.keys()):
        l_ism.append(
            [
                ism_pos[pos]["A"],
                ism_pos[pos]["T"],
                ism_pos[pos]["G"],
                ism_pos[pos]["C"],
            ]
        )
    ism = np.array(l_ism)
    ism_score = np.array([ism_pos[pos]["effect"] for pos in ism_pos])

    # All outputs need to be cropped to center +/- viewpoint_width/2
    return ism_pos, ism, ism_score


def get_saliency(model, seq):
    seq_tensor = tf.Variable(seq, dtype=tf.float32)
    with tf.GradientTape() as tape:
        tape.watch(seq_tensor)
        prediction = model(seq_tensor)
        output = tf.reduce_sum(prediction)

    gradients = tape.gradient(output, seq_tensor)[:, :, 1:5]
    gradient = gradients * seq_tensor[:, :, 1:5]
    gradient = tf.reduce_sum(gradient, axis=-1).numpy()[0]
    return gradient, gradients.numpy()[0]


def attribution(
    repo,
    model,
    track_name: str,
    chromosome: str,
    centerpoint: int,
    viewpoint_width: int = 0,
    batch: int = 0,
    skip_ism: bool = False,
) -> tuple:
    """
    Perform in-silico saturation mutagenesis (ISM) and calculate saliency
    of a given head over a given genomic region.

    Args:
        repo: Momics repository object
        model: Trained TensorFlow/Keras model for prediction
        track_name: Name of the track to analyze (e.g., "bw2")
        chromosome: Chromosome name (e.g., "chr1")
        centerpoint: Center point of the region of interest (int)
        viewpoint_width: Width of the region around the centerpoint to analyze (int)
        batch: Batch size for model prediction. If 0, no batching is applied (int)
        skip_ism: If True, skip ISM computation and return None for ISM-related outputs (bool)

    Returns:
        seq: One-hot encoded sequence of the region of interest
        data: Experimental data for the region of interest
        pred: Model predictions for the region of interest
        saliency: gradient x input of the model
        mutation_info: List of tuples (position, original_nuc, mutated_nuc)
        ism: a wx4 array, where w is the number of positions mutated, and the 4 columns are effects of mutating to A, T, G, C
        ism_score: the aggregated effect score for each position (average of the 3 mutations that are not the original nucleotide)
        ism_data: ism_score multiplied by the original experimental scores
    """

    ## Find coordinates
    in_width = model.input_shape[1]
    out_width = model.output_shape[1]
    s = centerpoint - in_width // 2
    e = centerpoint + in_width // 2
    if (e - s) > in_width:
        s = e - in_width
    elif (e - s) < in_width:
        e = s + in_width

    if viewpoint_width == 0:
        viewpoint_width = out_width
    os = centerpoint - viewpoint_width // 2
    oe = centerpoint + viewpoint_width // 2
    if (oe - os) > viewpoint_width:
        os = oe - viewpoint_width
    elif (oe - os) < viewpoint_width:
        oe = os + viewpoint_width

    coord_compute = f"{chromosome}:{s}-{e}"
    coord_out = f"{chromosome}:{os}-{oe}"

    ## Extract sequence and data for the ROI
    seq_compute = np.reshape(
        mutils.one_hot_encode(
            mmq.MomicsQuery(repo, coord_compute).query_sequence().seq["nucleotide"][coord_compute]  # type: ignore
        ),
        (1, in_width, 5),
    )
    data_compute = mmq.MomicsQuery(repo, coord_out).query_tracks(tracks=[track_name]).coverage[track_name][coord_out]  # type: ignore
    original_pred = model.predict(seq_compute, verbose=0)

    # Extract saliency
    saliency, _ = get_saliency(model, seq_compute)

    # Perform ISM
    if skip_ism:
        ism_pos, ism, ism_score, ism_data = None, None, None, None
    else:
        ism_pos, ism, ism_score = get_ISM(model, seq_compute, viewpoint_width, batch)
        ism_data = ism_score * data_compute

    return seq_compute, data_compute, original_pred, saliency, ism_pos, ism, ism_score, ism_data

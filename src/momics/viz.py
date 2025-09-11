from .query import MomicsQuery
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import Optional, Tuple


DEFAULT_ISM_CMAP = plt.get_cmap("afmhot_r")


def coverage(q: MomicsQuery, cmap: Optional[dict] = None, fig_size: Optional[Tuple[int, int]] = None) -> plt.Axes:
    """
    Plot the coverage of tracks queried over a single genomic window.

    Args:
        q (MomicsQuery): A MomicsQuery object containing the coverage data.
        cmap (dict, optional): A dictionary mapping track names to colors. Defaults to None.
        fig_size (Tuple[int, int], optional): Size of the figure. Defaults to (10, 1.5 * n_tracks).

    Returns:
        plt.Axes: The matplotlib Axes object containing the plot.

    Example:
        >>> repo = Momics("tests_data/S288c_MTL.momics")
        >>> q = MomicsQuery(repo, "IV:324469-332598").query_tracks(tracks=["MNASE", "H3K4ME3", "ATAC", "SUA7"])
        >>> fig = momics.viz.coverage(q)
    """
    if q.coverage is None:
        raise ValueError("Coverage data is not available in the MomicsQuery object.")

    cov = q.coverage
    gr = q.ranges
    tracks = cov.keys()
    if len(gr) != 1:
        raise ValueError("The query must contain exactly one genomic window.")

    chrom = gr.Chromosome[0]
    st = gr.Start[0]
    en = gr.End[0]
    idx = np.arange(st, en)
    cov_arr = np.stack([cov[track][f"{chrom}:{st}-{en}"] for track in tracks])

    n_tracks = len(tracks)
    if fig_size is None:
        fig_size = (10, int(1.5 * n_tracks))

    if cmap is not None:
        sns.set_palette(sns.color_palette([cmap[track] for track in tracks if track in cmap]))
    else:
        sns.set_palette(sns.color_palette("tab10", n_colors=len(tracks)))

    _, axes = plt.subplots(n_tracks, 1, figsize=fig_size, sharex=True)
    if n_tracks == 1:
        axes = [axes]
    for i, track in enumerate(tracks):
        ax = axes[i]
        coverage = cov_arr[i]
        sns.lineplot(x=idx, y=coverage, ax=ax, label=track, color=sns.color_palette()[i])
        ax.fill_between(idx, coverage, alpha=0.3, color=sns.color_palette()[i])
        ax.set_ylabel(track)
        ax.set_ylim(bottom=0)
        ax.set_xlim(st, en)
        ax.grid(True)
        if i == n_tracks - 1:
            ax.set_xlabel("Position")
        if i == 0:
            ax.set_title(f"{chrom}:{st}-{en}")

    plt.tight_layout()

    return axes[0].figure.axes


def aggrcoverage(
    q: MomicsQuery, ci: bool = True, cmap: Optional[dict] = None, fig_size: Optional[Tuple[int, int]] = None
) -> plt.Axes:
    """
    Plot the aggregate coverage of tracks queried over a number of genomic windows.

    Args:
        q (MomicsQuery): A MomicsQuery object containing the coverage data.
        ci (bool, optional): Whether to plot confidence intervals. Defaults to True.
        cmap (dict, optional): A dictionary mapping track names to colors. Defaults to None.
        fig_size (Tuple[int, int], optional): Size of the figure. Defaults to (10, 1.5 * n_tracks).

    Returns:
        plt.Axes: The matplotlib Axes object containing the plot.
    """
    if q.coverage is None:
        raise ValueError("Coverage data is not available in the MomicsQuery object.")

    tracks = q.coverage.keys()
    gr = q.ranges
    w = gr.End - gr.Start
    w = w.unique()
    if len(w) != 1:
        raise ValueError("The query must be performed over genomic regions of equal length.")

    w = w[0]
    st = -w // 2
    en = w // 2
    idx = np.arange(st, en)

    covs: dict = {x: {"mean": [], "ci": []} for x in tracks}
    for _, track_name in enumerate(tracks):
        scores = np.stack(list(q.coverage[track_name].values()))
        mean = scores.mean(axis=0)
        sd = scores.std(axis=0)
        covs[track_name]["mean"] = mean
        covs[track_name]["ci"] = 1.96 * sd / np.sqrt(scores.shape[0])

    if fig_size is None:
        fig_size = (7, 7)

    if cmap is not None:
        sns.set_palette(sns.color_palette([cmap[track] for track in tracks if track in cmap]))
    else:
        sns.set_palette(sns.color_palette("tab10", n_colors=len(tracks)))

    max_y = -99999

    plt.figure(figsize=fig_size)
    for i, track in enumerate(tracks):
        coverage = covs[track]["mean"]
        max_y = max(max_y, np.nanmax(coverage + covs[track]["ci"]))
        sns.lineplot(x=idx, y=coverage, label=track, color=sns.color_palette()[i])
        if ci:
            plt.fill_between(
                idx,
                coverage - covs[track]["ci"],
                coverage + covs[track]["ci"],
                alpha=0.3,
                color=sns.color_palette()[i],
                edgecolor=None,
            )
        plt.ylabel("Coverage")
        plt.ylim(bottom=0)
        plt.xlim(st, en)
        plt.grid(True)
        plt.xlabel("Position")

    plt.ylim(0, max_y * 1.1)

    plt.tight_layout()

    return plt.gca()


def heatcoverage(
    q: MomicsQuery, order_by: Optional[str] = None, cmap: Optional[dict] = None, fig_size: Optional[Tuple[int, int]] = None
) -> plt.Axes:
    """
    Plot heatmaps of queried coverage data over genomic regions.
    Each track is plotted in a separate heatmap.

    Args:
        q (MomicsQuery): A MomicsQuery object containing the coverage data.
        order_by (str, optional): Track name to order the regions by. Defaults to None.
        cmap (dict, optional): A dictionary mapping track names to colors. Defaults to None.
        fig_size (Tuple[int, int], optional): Size of the figure. Defaults to (10, 1.5 * n_tracks).

    Returns:
        plt.Axes: The matplotlib Axes object containing the plot.
    """
    if q.coverage is None:
        raise ValueError("Coverage data is not available in the MomicsQuery object.")

    tracks = q.coverage.keys()
    ntracks = len(tracks)
    gr = q.ranges
    if order_by is not None:
        if order_by not in tracks:
            raise ValueError(f"Track '{order_by}' not found in the queried tracks.")

    w = gr.End - gr.Start
    w = w.unique()
    if len(w) != 1:
        raise ValueError("The query must be performed over genomic regions of equal length.")

    w = w[0]
    st = -w // 2
    en = w // 2
    idx = np.arange(st, en)

    covs: dict = {x: None for x in tracks}
    for track_name in tracks:
        covs[track_name] = np.stack(list(q.coverage[track_name].values()))

    if order_by is not None:
        order_idx = np.argsort(covs[order_by].mean(axis=1))
        for track_name in tracks:
            covs[track_name] = covs[track_name][order_idx, :]

    if fig_size is None:
        fig_size = (7, 7)

    if cmap is not None:
        sns.set_palette(sns.color_palette([cmap[track] for track in tracks if track in cmap]))
    else:
        sns.set_palette(sns.color_palette("tab10", n_colors=ntracks))

    _, axes = plt.subplots(1, ntracks, figsize=fig_size)
    for i, track in enumerate(tracks):
        coverage = covs[track]
        cm = sns.light_palette(sns.color_palette()[i], as_cmap=True)
        sns.heatmap(coverage, robust=True, cbar=True, cmap=cm, ax=axes[i])
        axes[i].set_ylabel("")
        axes[i].set_xlabel("Position")
        axes[i].set_ylim(bottom=0, top=coverage.shape[0])
        axes[i].set_yticks([])
        axes[i].set_xticks([0, len(idx) // 2, len(idx) - 1])
        axes[i].set_xticklabels([idx[0] + 1, 0, idx[len(idx) - 1] + 1], rotation=0)
        axes[i].grid(True)
        axes[i].set_title(track)
        cbar = axes[i].collections[0].colorbar
        cbar.outline.set_edgecolor("black")
        cbar.outline.set_linewidth(0.5)
        for _, spine in axes[i].spines.items():
            spine.set_visible(True)
            spine.set_linewidth(0.5)
            spine.set_color("black")

    plt.tight_layout()

    return plt.gca()


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

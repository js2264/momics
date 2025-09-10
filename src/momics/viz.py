from .query import MomicsQuery
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import Optional, Tuple


def coverage(q: MomicsQuery, fig_size: Optional[Tuple[int, int]] = None) -> plt.Axes:
    """
    Plot the coverage of tracks queried over a single genomic window.

    Args:
        q (MomicsQuery): A MomicsQuery object containing the coverage data.
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

    _, axes = plt.subplots(n_tracks, 1, figsize=fig_size, sharex=True)
    if n_tracks == 1:
        axes = [axes]
    for i, track in enumerate(tracks):
        ax = axes[i]
        coverage = cov_arr[i]
        sns.lineplot(x=idx, y=coverage, ax=ax, label=track, color=sns.color_palette("tab10")[i])
        ax.fill_between(idx, coverage, alpha=0.3, color=sns.color_palette("tab10")[i])
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

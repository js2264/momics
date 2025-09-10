import matplotlib
import matplotlib.pyplot as plt
import pytest
from momics.momics import Momics
from momics.query import MomicsQuery
from momics import utils as mutils
from momics import viz

matplotlib.use("Agg")


@pytest.mark.order(10)
def test_coverage_plot():
    """Test the coverage plotting function for a single genomic window."""
    mom = Momics("tests_data/test.momics")
    grs = mutils.parse_ucsc_coordinates(["I:990-1010", "II:200-300"])
    q0 = MomicsQuery(mom, "I:990-1010").query_tracks(tracks=["ATAC", "SCC1"])
    q1 = MomicsQuery(mom, grs).query_tracks(tracks=["ATAC", "SCC1"])
    axes = viz.coverage(q0)
    assert axes is not None
    assert len(axes) == 2
    for ax in axes:
        assert ax.get_xlabel() == "Position" or ax.get_xlabel() == ""
        assert ax.get_xlim() == (990, 1010)
        assert ax.get_ylim()[0] == 0

    axes_custom = viz.coverage(q0, fig_size=(15, 5))
    assert axes_custom is not None

    with pytest.raises(ValueError, match="must contain exactly one genomic window"):
        viz.coverage(q1)

    plt.close("all")


@pytest.mark.order(10)
def test_coverage_plot_single_track():
    """Test coverage plotting with a single track."""
    mom = Momics("tests_data/test.momics")
    q = MomicsQuery(mom, "I:990-1010").query_tracks(tracks=["ATAC"])
    axes = viz.coverage(q)
    assert axes is not None
    assert len(axes) == 1
    ax = axes[0]
    assert ax.get_xlabel() == "Position"
    assert ax.get_xlim() == (990, 1010)
    assert ax.get_ylim()[0] == 0
    plt.close("all")


@pytest.mark.order(10)
def test_coverage_plot_errors():
    """Test error handling in coverage plotting."""
    mom = Momics("tests_data/test.momics")
    q = MomicsQuery(mom, "I:990-1010")
    with pytest.raises(ValueError, match="Coverage data is not available"):
        viz.coverage(q)


@pytest.mark.order(10)
def test_aggrcoverage_plot() -> None:
    """Test the aggregate coverage plotting function."""
    mom = Momics("tests_data/test.momics")
    ranges = mutils.parse_ucsc_coordinates(["I:0-10", "I:20-30"])
    q = MomicsQuery(mom, ranges).query_tracks(tracks=["ATAC", "SCC1"])
    ax = viz.aggrcoverage(q)
    assert ax is not None
    assert ax.get_xlabel() == "Position"
    assert ax.get_ylabel() == "Coverage"
    assert ax.get_xlim() == (-5, 5)
    assert ax.get_ylim()[0] == 0
    plt.close("all")


@pytest.mark.order(10)
def test_aggrcoverage_plot_no_ci():
    """Test aggregate coverage plotting without confidence intervals."""
    mom = Momics("tests_data/test.momics")
    ranges = mutils.parse_ucsc_coordinates(["I:0-10", "I:20-30"])
    q = MomicsQuery(mom, ranges).query_tracks(tracks=["ATAC", "SCC1"])
    ax = viz.aggrcoverage(q, ci=False)
    assert ax is not None
    plt.close("all")


@pytest.mark.order(10)
def test_aggrcoverage_plot_custom_cmap():
    """Test aggregate coverage plotting with custom colormap."""
    mom = Momics("tests_data/test.momics")
    ranges = mutils.parse_ucsc_coordinates(["I:0-10", "I:20-30"])
    q = MomicsQuery(mom, ranges).query_tracks(tracks=["ATAC", "SCC1"])
    custom_cmap = {"ATAC": "red", "SCC1": "blue"}
    ax = viz.aggrcoverage(q, cmap=custom_cmap)
    assert ax is not None
    ax_custom = viz.aggrcoverage(q, fig_size=(8, 6))
    assert ax_custom is not None
    plt.close("all")


@pytest.mark.order(10)
def test_aggrcoverage_plot_errors():
    """Test error handling in aggregate coverage plotting."""
    mom = Momics("tests_data/test.momics")
    ranges = mutils.parse_ucsc_coordinates(["I:0-10", "I:20-30"])
    q = MomicsQuery(mom, ranges)
    with pytest.raises(ValueError, match="Coverage data is not available"):
        viz.aggrcoverage(q)

    ranges_unequal = mutils.parse_ucsc_coordinates(["I:0-10", "I:20-35"])
    q_unequal = MomicsQuery(mom, ranges_unequal).query_tracks(tracks=["ATAC"])
    with pytest.raises(ValueError, match="equal length"):
        viz.aggrcoverage(q_unequal)


@pytest.mark.order(10)
def test_viz_functions_return_types():
    """Test that visualization functions return expected types."""
    mom = Momics("tests_data/test.momics")
    q_single = MomicsQuery(mom, "I:990-1010").query_tracks(tracks=["ATAC"])
    result_coverage = viz.coverage(q_single)
    assert isinstance(result_coverage, list)
    ranges = mutils.parse_ucsc_coordinates(["I:0-10", "I:20-30"])
    q_multi = MomicsQuery(mom, ranges).query_tracks(tracks=["ATAC"])
    result_aggr = viz.aggrcoverage(q_multi)
    assert hasattr(result_aggr, "get_xlabel")
    plt.close("all")


@pytest.mark.order(10)
def test_viz_functions_with_empty_coverage():
    """Test visualization functions behavior with tracks containing zeros."""
    mom = Momics("tests_data/test.momics")
    q = MomicsQuery(mom, "I:5000-5020").query_tracks(tracks=["ATAC", "SCC1"])
    axes = viz.coverage(q)
    assert axes is not None
    ranges = mutils.parse_ucsc_coordinates(["I:5000-5020", "I:6000-6020"])
    q_multi = MomicsQuery(mom, ranges).query_tracks(tracks=["ATAC", "SCC1"])
    ax = viz.aggrcoverage(q_multi)
    assert ax is not None
    plt.close("all")


@pytest.mark.order(10)
def test_heatcoverage():
    """Test heatcoverage plotting function."""
    mom = Momics("tests_data/test.momics")
    bins = mom.bins(1000, 1000, cut_last_bin_out=True)["I"]
    q = MomicsQuery(mom, bins).query_tracks(tracks=["ATAC", "SCC1"])
    axes = viz.heatcoverage(q, order_by="SCC1", cmap={"ATAC": "red", "SCC1": "blue"})
    assert axes is not None
    plt.close("all")

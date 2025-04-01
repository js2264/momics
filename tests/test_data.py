from typing import Generator, Iterator
import pytest

from momics import momics
from momics.streamer import MomicsStreamer
from momics.dataset import MomicsDataset


@pytest.mark.order(99)
def test_streamer(momics_path: str):
    mom = momics.Momics(momics_path)

    b = mom.bins(10, 21, cut_last_bin_out=True)
    with pytest.raises(ValueError, match=r".*not found in momics repository."):
        MomicsStreamer(mom, b, features=["CH1", "bw2"])

    rg = MomicsStreamer(mom, b, features="bw2", batch_size=1000, silent=False)
    assert isinstance(rg.generator(), Generator)
    assert isinstance(iter(rg), Iterator)

    n = next(iter(rg))
    assert len(n) == 1
    assert n[0].shape == (1000, 10, 1)
    n = next(rg)
    assert len(n) == 1
    assert n[0].shape == (1000, 10, 1)

    for n in iter(rg):
        print(rg.batch_index, " / ", rg.num_batches)
        assert isinstance(n, tuple)
    with pytest.raises(StopIteration):
        next(rg)
    assert rg.batch_index == rg.num_batches

    rg = MomicsStreamer(mom, b, features=["bw3", "bw2"], batch_size=1000)
    n = next(rg)
    assert len(n) == 2
    assert n[0].shape == (1000, 10, 1)
    assert n[1].shape == (1000, 10, 1)
    n = next(rg)
    assert len(n) == 2
    assert n[0].shape == (1000, 10, 1)
    assert n[1].shape == (1000, 10, 1)

    rg = MomicsStreamer(mom, b, features=["nucleotide", "bw2"], batch_size=10)
    n = next(rg)
    assert len(n) == 2
    assert n[0].shape == (10, 10, 4)
    assert n[1].shape == (10, 10, 1)
    n = next(rg)
    assert len(n) == 2
    assert n[0].shape == (10, 10, 4)
    assert n[1].shape == (10, 10, 1)


@pytest.mark.order(99)
def test_dataset(momics_path: str):
    mom = momics.Momics(momics_path)
    b = mom.bins(10, 21)

    with pytest.raises(ValueError, match="All ranges must have the same width"):
        MomicsDataset(mom, b, "CH0", "CH1")

    b = mom.bins(10, 21, cut_last_bin_out=True)
    with pytest.raises(ValueError, match=r"Target size must be smaller.*"):
        MomicsDataset(mom, b, "bw3", "bw2", target_size=1000000)

    rg = MomicsDataset(mom, b, "bw3", "bw2", target_size=2, batch_size=10)
    n = next(iter(rg))
    assert n[0][0].shape == (10, 10, 1)
    assert n[1][0].shape == (10, 2, 1)

    rg = MomicsDataset(mom, b, "nucleotide", "bw2", target_size=2, batch_size=10)
    n = next(iter(rg))
    assert n[0][0].shape == (10, 10, 4)
    assert n[1][0].shape == (10, 2, 1)
    with pytest.raises(IndexError, match=r"tuple index out of range"):
        _ = n[0][1].shape
    with pytest.raises(IndexError, match=r"tuple index out of range"):
        _ = n[1][1].shape

    # Only features, list of features
    rg = MomicsDataset(mom, b, ["nucleotide", "bw2"], target_size=2, batch_size=10)
    n = next(iter(rg))
    assert n[0][0].shape == (10, 10, 4)
    assert n[1][0].shape == (10, 10, 1)
    with pytest.raises(IndexError, match=r"tuple index out of range"):
        _ = n[2][0].shape

    # List of features and single target
    rg = MomicsDataset(mom, b, ["nucleotide", "bw2"], target="bw3", target_size=2, batch_size=10)
    n = next(iter(rg))
    assert n[0][0][0].shape == (10, 10, 4)
    assert n[0][1][0].shape == (10, 10, 1)
    assert n[1][0].shape == (10, 2, 1)
    with pytest.raises(IndexError, match=r"tuple index out of range"):
        _ = n[1][1].shape

    # List of features and list of targets
    rg = MomicsDataset(mom, b, ["nucleotide", "bw2"], target=["bw3", "bw2"], target_size=2, batch_size=10)
    n = next(iter(rg))
    assert n[0][0][0].shape == (10, 10, 4)
    assert n[0][1][0].shape == (10, 10, 1)
    assert n[1][0][0].shape == (10, 2, 1)
    assert n[1][1][0].shape == (10, 2, 1)

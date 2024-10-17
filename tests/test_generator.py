import pytest

import momics
from momics.dataloader import RangeDataLoader


@pytest.mark.order(99)
def test_generator_init(momics_path: str):
    mom = momics.Momics(momics_path)
    b = mom.bins(10, 21)

    with pytest.raises(ValueError, match="All ranges must have the same width"):
        RangeDataLoader(mom, b, "CH0", "CH1")

    b = mom.bins(10, 20, cut_last_bin_out=True)
    with pytest.raises(ValueError, match=r".*not found in momics repository"):
        RangeDataLoader(mom, b, "CH0", "CH1")
    with pytest.raises(ValueError, match=r".*not found in momics repository"):
        RangeDataLoader(mom, b, "bw3", "CH1")
    with pytest.raises(ValueError, match=r"Target size must be smaller than the features width"):
        RangeDataLoader(mom, b, "bw3", "bw2", target_size=1000000)

    rg = RangeDataLoader(mom, b, "bw3", "bw2", target_size=2, batch_size=1000)
    n = rg[0]
    assert n[0].shape == (1000, 10, 1)
    assert n[1].shape == (1000, 2, 1)

    rg = RangeDataLoader(mom, b, "nucleotide", "bw2", target_size=2, batch_size=1000)
    n = rg[0]
    assert n[0].shape == (1000, 10, 4)
    assert n[1].shape == (1000, 2, 1)

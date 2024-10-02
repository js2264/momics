import numpy as np
import pandas as pd
import pytest
import tiledb
from pybedtools import BedTool

import momics
from momics import utils
from momics.multirangequery import MultiRangeQuery


@pytest.mark.order(1)
def test_Momics_init(momics_path: str):

    x = momics.Momics(momics_path)
    assert x.path == momics_path


@pytest.mark.order(1)
def test_Momics_add_genome(momics_path: str, bw1: str):
    mom = momics.Momics(momics_path)

    assert mom.chroms().empty

    with pytest.raises(ValueError, match=r"Please fill out `chroms` table first."):
        mom.add_tracks({"bw1": bw1})

    chroms = utils.get_chr_lengths(bw1)
    mom.add_chroms(chroms)
    out = pd.DataFrame(
        {
            "chrom_index": [0, 1, 2],
            "chrom": ["I", "II", "III"],
            "length": [10000, 20000, 30000],
        }
    )
    assert mom.chroms().__eq__(out).all().all()
    with pytest.raises(
        ValueError, match=r"`chroms` table has already been filled out."
    ):
        mom.add_chroms(chroms)


@pytest.mark.order(1)
def test_Momics_add_tracks(momics_path: str, bw1: str, bw2: str):
    mom = momics.Momics(momics_path)

    assert mom.tracks().empty
    mom.add_tracks({"bw1": bw1})
    out = pd.DataFrame(
        {
            "idx": [0],
            "label": ["bw1"],
            "path": [bw1],
        }
    )
    assert mom.tracks().__eq__(out).all().all()
    with pytest.raises(ValueError, match=r".*already present in `tracks` table"):
        mom.add_tracks({"bw1": bw1})
    with pytest.raises(Exception, match=r".*do not have identical chromomosome.*"):
        mom.add_tracks({"bw2": bw2})
    mom.add_tracks({"bw2": bw1})
    out = pd.DataFrame(
        {
            "idx": [0, 1],
            "label": ["bw1", "bw2"],
            "path": [bw1, bw1],
        }
    )
    assert mom.tracks().__eq__(out).all().all()
    print(mom.tracks())


@pytest.mark.order(1)
def test_Momics_add_track(momics_path: str, bw1: str, bw2: str):
    mom = momics.Momics(momics_path)
    chroms = mom.chroms()
    coverage = {
        chrom: np.random.rand(length) for i, (idx, chrom, length) in chroms.iterrows()
    }
    with pytest.raises(ValueError, match=r".*already present in `tracks` table"):
        mom.add_track(coverage, "bw1")

    mom.add_track(coverage, "custom")
    print(mom.tracks())
    out = pd.DataFrame(
        {
            "idx": [0, 1, 2],
            "label": ["bw1", "bw2", "custom"],
            "path": [bw1, bw1, "custom"],
        }
    )
    assert mom.tracks().__eq__(out).all().all()
    print(mom.tracks())


@pytest.mark.order(1)
def test_Momics_add_seq(momics_path: str, fa1: str, fa2: str):
    mom = momics.Momics(momics_path)

    with pytest.raises(Exception, match=r".*do not have identical chromomosome.*"):
        mom.add_sequence(fa2)

    mom.add_sequence(fa1)

    with pytest.raises(tiledb.cc.TileDBError, match=r"already exists"):
        mom.add_sequence(fa2)

    print(mom.seq())


@pytest.mark.order(1)
def test_Momics_remove_tracks(momics_path: str, bw1: str, bw2: str, bed1: str):
    mom = momics.Momics(momics_path)
    mom.add_tracks({"bw3": bw1})
    mom.add_tracks({"bw4": bw1})
    mom.remove_track("bw1")
    print(mom.tracks())
    out = pd.DataFrame(
        {
            "idx": [0, 1, 2, 3, 4],
            "label": ["None", "bw2", "custom", "bw3", "bw4"],
            "path": ["None", bw1, "custom", bw1, bw1],
        }
    )
    print(out)
    assert mom.tracks().__eq__(out).all().all()
    q = MultiRangeQuery(mom, "I:991-1010").query_tracks()
    assert list(q.coverage.keys()) == ["bw2", "custom", "bw3", "bw4"]
    bed = BedTool(bed1).to_dataframe()
    q = MultiRangeQuery(mom, bed).query_tracks()
    assert list(q.coverage.keys()) == ["bw2", "custom", "bw3", "bw4"]


@pytest.mark.order(2)
def test_Momics_binnify(momics_path: str):
    mom = momics.Momics(momics_path)
    q = mom.bins(width=1000, step=1000)
    assert q.shape == (60, 3)

import numpy as np
import pandas as pd
import pytest

import momics
from momics import utils
from momics.multirangequery import MultiRangeQuery


def test_Momics_init(momics_path: str):
    with pytest.raises(OSError, match=r"Momics repository not found."):
        momics.Momics("asvasdvasdv", create=False)

    momics.Momics(momics_path)
    x = momics.Momics(momics_path, create=False)

    assert x.path == momics_path

    with pytest.raises(OSError, match=r".* repository already exist found"):
        momics.Momics(momics_path, create=True)


def test_Momics_add_genome(momics_path: str, bw1: str):
    mom = momics.Momics(momics_path, create=False)

    assert mom.chroms().empty

    with pytest.raises(ValueError, match=r"Please fill out `chroms` table first."):
        mom.add_tracks({"bw1": bw1})

    chroms = utils.get_chr_lengths(bw1)
    mom.add_chroms(chroms)
    out = pd.DataFrame(
        {
            "chrom_index": [0, 1, 2],
            "chr": ["I", "II", "III"],
            "length": [10000, 20000, 30000],
        }
    )
    assert mom.chroms().__eq__(out).all().all()
    with pytest.raises(
        ValueError, match=r"`chroms` table has already been filled out."
    ):
        mom.add_chroms(chroms)


def test_Momics_add_tracks(momics_path: str, bw1: str, bw2: str):
    mom = momics.Momics(momics_path, create=False)

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


def test_Momics_add_seq(momics_path: str, fa1: str, fa2: str):
    mom = momics.Momics(momics_path, create=False)

    with pytest.raises(Exception, match=r".*do not have identical chromomosome.*"):
        mom.add_sequence(fa2)

    mom.add_sequence(fa1)

    with pytest.raises(ValueError, match=r"Sequence already added to the repository"):
        mom.add_sequence(fa2)

    print(mom.seq())


def test_Momics_query_tracks(momics_path: str, bed1: str):
    mom = momics.Momics(momics_path, create=False)
    q = MultiRangeQuery(mom, "I:991-1010").query_tracks()
    assert len(q.coverage) == 2
    assert len(q.coverage["bw1"]["I:991-1010"]) == 20

    q = MultiRangeQuery(mom, "I").query_tracks()
    assert len(q.coverage) == 2
    assert len(q.coverage["bw1"]["I:1-10000"]) == 10000

    q = MultiRangeQuery(mom, "I").query_tracks()
    assert len(q.coverage) == 2
    assert len(q.coverage["bw1"]["I:1-10000"]) == 10000

    # q = MultiRangeQuery(mom, "I").query_tracks(with_seq=True)
    # assert len(q.coverage) == 3

    bed = utils.import_bed_file(bed1)
    q = MultiRangeQuery(mom, bed).query_tracks()
    assert q.coverage.keys().__eq__(["bw2", "bw3", "bw4"])


def test_Momics_query_seqs(momics_path: str):
    mom = momics.Momics(momics_path, create=False)
    q = mom.query_sequence("I:991-1010")
    assert len(q) == 20

    q = mom.query_sequence("I")
    assert len(q) == 10000


def test_Momics_remove_tracks(momics_path: str, bw1: str, bw2: str, bed1: str):
    mom = momics.Momics(momics_path, create=False)
    mom.add_tracks({"bw3": bw1})
    mom.add_tracks({"bw4": bw1})
    mom.remove_track("bw1")
    out = pd.DataFrame(
        {
            "idx": [0, 1, 2, 3],
            "label": ["None", "bw2", "bw3", "bw4"],
            "path": ["None", bw1, bw1, bw1],
        }
    )
    assert mom.tracks().__eq__(out).all().all()
    q = MultiRangeQuery(mom, "I:991-1010").query_tracks()
    assert q.coverage.keys().__eq__(["bw2", "bw3", "bw4"])
    bed = utils.import_bed_file(bed1)
    q = MultiRangeQuery(mom, bed).query_tracks()
    assert q.coverage.keys().__eq__(["bw2", "bw3", "bw4"])


def test_Momics_binnify(momics_path: str):
    mom = momics.Momics(momics_path, create=False)
    q = mom.bins(width=1000, step=1000)
    assert q.shape == (60, 3)

import pandas as pd
from pybedtools import BedTool
import pytest

import momics
from momics.multirangequery import MultiRangeQuery


@pytest.mark.order(2)
def test_multirangequery_tracks(momics_path: str, bed1: str):
    mom = momics.Momics(momics_path, create=False)
    q = MultiRangeQuery(mom, "I:991-1010").query_tracks()
    assert len(q.coverage) == 3
    assert len(q.coverage["bw2"]["I:991-1010"]) == 20
    assert q.to_df()["chrom"].__eq__(pd.Series(["I"] * 20)).all()

    q = MultiRangeQuery(mom, "I").query_tracks()
    assert len(q.coverage) == 3
    assert len(q.coverage["bw2"]["I:1-10000"]) == 10000

    q = MultiRangeQuery(mom, "I").query_tracks()
    assert len(q.coverage) == 3
    assert len(q.coverage["bw2"]["I:1-10000"]) == 10000

    bed = BedTool(bed1).to_dataframe()
    q = MultiRangeQuery(mom, bed).query_tracks()
    assert q.coverage.keys().__eq__(["bw2", "bw3", "bw4"])

    q = MultiRangeQuery(mom, bed).query_tracks(threads=4)
    assert q.coverage.keys().__eq__(["bw2", "bw3", "bw4"])


@pytest.mark.order(2)
def test_multirangequery_seq(momics_path: str, bed1: str):
    mom = momics.Momics(momics_path, create=False)
    q = MultiRangeQuery(mom, "I:1-10").query_sequence()
    assert len(q.seq) == 1
    assert q.seq["seq"]["I:1-10"] == "ATCGATCGAT"

    bed = BedTool(bed1).to_dataframe()
    q = MultiRangeQuery(mom, bed).query_sequence()
    assert len(q.seq["seq"]) == 3
    assert q.seq["seq"]["I:1-10"] == "ATCGATCGAT"

    assert q.to_fa()[0].id == "I:1-10"

    q = MultiRangeQuery(mom, "I").query_sequence()
    assert len(q.seq["seq"]["I:1-10000"]) == 10000

    q = MultiRangeQuery(mom, "I").query_sequence(4)
    assert len(q.seq["seq"]["I:1-10000"]) == 10000

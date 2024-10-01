import json
import os
from pathlib import Path
import pickle
import tempfile
import numpy as np
import pandas as pd
from pybedtools import BedTool
import pytest

import momics
from momics.multirangequery import MultiRangeQuery


@pytest.mark.order(2)
def test_multirangequery_tracks(momics_path: str, bed1: str):
    mom = momics.Momics(momics_path, create=False)
    q = MultiRangeQuery(mom, "I:991-1010").query_tracks()
    assert len(q.coverage) == 4
    assert len(q.coverage["bw2"]["I:991-1010"]) == 20
    assert q.to_df()["chrom"].__eq__(pd.Series(["I"] * 20)).all()

    q = MultiRangeQuery(mom, "I").query_tracks()
    assert len(q.coverage) == 4
    assert len(q.coverage["bw2"]["I:1-10000"]) == 10000

    q = MultiRangeQuery(mom, "I").query_tracks()
    assert len(q.coverage) == 4
    assert len(q.coverage["bw2"]["I:1-10000"]) == 10000

    bed = BedTool(bed1).to_dataframe()
    q = MultiRangeQuery(mom, bed).query_tracks()
    assert q.coverage.keys().__eq__(["bw2", "custom", "bw3", "bw4"])

    q = MultiRangeQuery(mom, bed).query_tracks(threads=4)
    assert q.coverage.keys().__eq__(["bw2", "custom", "bw3", "bw4"])


@pytest.mark.order(2)
def test_multirangequery_seq(momics_path: str):
    mom = momics.Momics(momics_path, create=False)
    q = MultiRangeQuery(mom, "I:1-10").query_sequence()
    assert len(q.seq) == 1
    assert q.seq["seq"]["I:1-10"] == "ATCGATCGAT"


@pytest.mark.order(2)
def test_multirangequery_seq2(momics_path: str, bed1: str):
    mom = momics.Momics(momics_path, create=False)
    q = MultiRangeQuery(mom, "I:1-10").query_sequence()
    bed = BedTool(bed1).to_dataframe()
    q = MultiRangeQuery(mom, bed).query_sequence()
    assert len(q.seq["seq"]) == 3
    assert q.seq["seq"]["I:1-10"] == "ATCGATCGAT"
    assert q.to_fa()[0].id == "I:1-10"


@pytest.mark.order(2)
def test_multirangequery_seq3(momics_path: str):
    mom = momics.Momics(momics_path, create=False)
    q = MultiRangeQuery(mom, "I:1-10").query_sequence()
    q = MultiRangeQuery(mom, "I").query_sequence()
    assert len(q.seq["seq"]["I:1-10000"]) == 10000


@pytest.mark.order(2)
def test_multirangequery_seq4(momics_path: str, bed1: str):
    mom = momics.Momics(momics_path, create=False)
    q = MultiRangeQuery(mom, "I:1-10").query_sequence()
    bed = BedTool(bed1).to_dataframe()
    q = MultiRangeQuery(mom, bed).query_sequence()
    assert q.seq["seq"]["I:1-10"] == "ATCGATCGAT"


@pytest.mark.order(2)
def test_multirangequery_seq5(momics_path: str, bed1: str):
    mom = momics.Momics(momics_path, create=False)
    bed = BedTool(bed1).to_dataframe()
    q = MultiRangeQuery(mom, "I:1-10").query_sequence()
    q = MultiRangeQuery(mom, bed).query_sequence(threads=4)
    assert q.seq["seq"]["I:1-10"] == "ATCGATCGAT"


@pytest.mark.order(2)
@pytest.fixture
def temp_json_file():
    """Fixture to create a temporary JSON for testing."""
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_json_file = Path(temp_file.name)
        yield temp_json_file
    temp_json_file.unlink()


@pytest.mark.order(2)
@pytest.fixture
def temp_npz_file():
    """Fixture to create a temporary NPZ for testing."""
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_npz_file = Path(temp_file.name)
        yield temp_npz_file
    temp_npz_file.unlink()


def test_to_json_npz(momics_path: str, temp_json_file: Path, temp_npz_file: Path):
    mom = momics.Momics(momics_path, create=False)
    q = MultiRangeQuery(mom, "I:1-10").query_sequence().query_tracks()

    q.to_json(temp_json_file)
    q.to_npz(temp_npz_file)

    with open(temp_json_file, "r") as json_file:
        data = json.load(json_file)
    assert data.keys().__eq__(["bw2", "custom", "bw3", "bw4", "seq"])
    assert data["seq"]["I:1-10"] == "ATCGATCGAT"

    data = np.load(temp_npz_file, allow_pickle=True)
    seq = pickle.loads(data["seq"])
    cov = pickle.loads(data["coverage"])
    assert seq["seq"]["I:1-10"] == "ATCGATCGAT"
    assert cov.keys().__eq__(["bw2", "custom", "bw3", "bw4"])
    assert cov["bw2"].keys().__eq__(["I:1-10"])

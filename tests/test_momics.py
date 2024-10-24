import collections
from pathlib import Path
import numpy as np
import pyranges as pr
import pandas as pd
import pytest
import tiledb

from momics import momics
from momics import utils
from momics.query import MomicsQuery


@pytest.mark.order(1)
def test_Momics_init(momics_path: str):

    x = momics.Momics(momics_path)
    assert x.path == momics_path


@pytest.mark.order(1)
def test_Momics_ingest_genome(momics_path: str, bw1: str):
    mom = momics.Momics(momics_path)

    assert mom.chroms().empty

    with pytest.raises(ValueError, match=r"Please fill out `chroms` table first."):
        mom.ingest_tracks({"bw1": bw1})

    chroms = utils.get_chr_lengths(bw1)
    mom.ingest_chroms(chroms)
    out = pd.DataFrame(
        {
            "chrom_index": [0, 1, 2],
            "chrom": ["I", "II", "III"],
            "length": [10000, 20000, 30000],
        }
    )
    assert mom.chroms().__eq__(out).all().all()
    with pytest.raises(ValueError, match=r"`chroms` table has already been filled out."):
        mom.ingest_chroms(chroms)


@pytest.mark.order(1)
def test_Momics_ingest_tracks(momics_path: str, bw1: str, bw2: str):
    mom = momics.Momics(momics_path)

    assert mom.tracks().empty
    mom.ingest_tracks({"bw1": bw1}, tile=10000)
    out = pd.DataFrame(
        {
            "idx": [0],
            "label": ["bw1"],
            "path": [bw1],
        }
    )
    assert mom.tracks().__eq__(out).all().all()
    with pytest.raises(ValueError, match=r".*already present in `tracks` table"):
        mom.ingest_tracks({"bw1": bw1})
    with pytest.raises(Exception, match=r".*do not have identical chromomosome.*"):
        mom.ingest_tracks({"bw2": bw2})
    mom.ingest_tracks({"bw2": bw1})
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
def test_Momics_alternative_ingest_tracks(momics_path2: str, bw1: str, bw2: str):
    mom = momics.Momics(momics_path2)
    chroms = utils.get_chr_lengths(bw1)
    momics.BULK_OVERWRITE = False
    mom.ingest_chroms(chroms)
    mom.ingest_tracks({"bw1": bw1})
    mom.ingest_tracks({"bw2": bw1})
    out = pd.DataFrame(
        {
            "idx": [0, 1],
            "label": ["bw1", "bw2"],
            "path": [bw1, bw1],
        }
    )
    assert mom.tracks().__eq__(out).all().all()
    print(mom.tracks())
    momics.BULK_OVERWRITE = True


@pytest.mark.order(1)
def test_Momics_ingest_track(momics_path: str, bw1: str, bw2: str):
    mom = momics.Momics(momics_path)
    chroms = mom.chroms()
    coverage = {chrom: np.random.rand(length) for i, (idx, chrom, length) in chroms.iterrows()}
    with pytest.raises(ValueError, match=r".*already present in `tracks` table"):
        mom.ingest_track(coverage, "bw1")

    mom.ingest_track(coverage, "custom")
    print(mom.tracks())
    out = pd.DataFrame(
        {
            "idx": [0, 1, 2],
            "label": ["bw1", "bw2", "custom"],
            # "path": [bw1, bw1, "custom"],
        }
    )
    assert mom.tracks().iloc[:, 0:2].__eq__(out).all().all()
    print(mom.tracks())


@pytest.mark.order(1)
def test_Momics_recover_track(momics_path: str):
    mom = momics.Momics(momics_path)
    print(mom.path)

    with pytest.raises(ValueError, match=r".*not found"):
        mom.tracks("bw1323")

    cov = mom.tracks("bw2")

    chrom_sizes = {"I": 10000, "II": 20000, "III": 30000}
    act = {chrom: [0] * length for chrom, length in chrom_sizes.items()}
    for chrom, size in chrom_sizes.items():
        intervals = [(i, i + 1000, i / 100000) for i in range(0, size, 1000)]
        x = [[v] * 1000 for (_, _, v) in intervals]
        arr = np.array([item for sublist in x for item in sublist], dtype=np.float32)
        act[chrom] = arr  # type: ignore

    for chrom in chrom_sizes.keys():
        assert cov[chrom].__eq__(act[chrom]).all()


@pytest.mark.order(1)
def test_Momics_ingest_seq(momics_path: str, fa1: str, fa2: str):
    mom = momics.Momics(momics_path)

    with pytest.raises(Exception, match=r".*do not have identical chromomosome.*"):
        mom.ingest_sequence(fa2)

    mom.ingest_sequence(fa1, tile=10000)

    with pytest.raises(tiledb.cc.TileDBError, match=r"already exists"):
        mom.ingest_sequence(fa2)

    print(mom.seq())

    assert mom.seq().shape == (3, 4)

    with pytest.raises(ValueError, match=r"Selected attribute does not exist.*"):
        mom.seq("csadc")
    assert isinstance(mom.seq("I"), str)
    assert len(mom.seq("II")) == 20000


@pytest.mark.order(1)
def test_Momics_remove_tracks(momics_path: str, bw1: str, bw2: str, bed1: str):
    mom = momics.Momics(momics_path)
    mom.ingest_tracks({"bw3": bw1})
    mom.ingest_tracks({"bw4": bw1})
    mom.remove_track("bw1")
    print(mom.tracks())
    out = pd.DataFrame(
        {
            "idx": [0, 1, 2, 3, 4],
            "label": ["None", "bw2", "custom", "bw3", "bw4"],
            # "path": ["None", bw1, "custom", bw1, bw1],
        }
    )
    print(out)
    assert mom.tracks().iloc[:, 0:2].__eq__(out).all().all()
    q = MomicsQuery(mom, "I:991-1010").query_tracks()
    assert list(q.coverage.keys()) == ["bw2", "custom", "bw3", "bw4"]
    bed = pr.read_bed(bed1)
    q = MomicsQuery(mom, bed).query_tracks()
    assert list(q.coverage.keys()) == ["bw2", "custom", "bw3", "bw4"]


@pytest.mark.order(2)
def test_Momics_binnify(momics_path: str):
    mom = momics.Momics(momics_path)
    q = mom.bins(width=1000, stride=1000)
    assert q.df.shape == (60, 3)


@pytest.mark.order(2)
def test_Momics_consolidate(momics_path: str):
    mom = momics.Momics(momics_path)
    s = mom.size()
    x = mom.consolidate(vacuum=True)
    assert x
    assert mom.size() < s


@pytest.mark.order(2)
def test_Momics_manifest(momics_path: str):
    mom = momics.Momics(momics_path)
    man = mom.manifest()
    assert isinstance(man, collections.defaultdict)


@pytest.mark.order(2)
def test_Momics_features(momics_path: str):
    mom = momics.Momics(momics_path)
    assert mom.features().empty

    sets = {
        "ft1": mom.bins(1000, 2000, cut_last_bin_out=True),
        "ft2": mom.bins(2, 24, cut_last_bin_out=True),
    }
    mom.ingest_features(sets, tile=10000)
    out = pd.DataFrame(
        {
            "idx": [0, 1],
            "label": ["ft1", "ft2"],
            "n": [30, 2501],
        }
    )
    assert mom.features().__eq__(out).all().all()

    with pytest.raises(ValueError, match=r".*already present in `features` table"):
        mom.ingest_features({"ft1": mom.bins(1000, 2000, cut_last_bin_out=True)}, tile=10000)

    sets = {
        "ft3": mom.bins(1000, 2000, cut_last_bin_out=True),
        "ft4": mom.bins(2, 24, cut_last_bin_out=True),
    }
    mom.ingest_features(sets)
    assert mom.features().shape == (4, 3)

    ft1 = mom.bins(1000, 2000, cut_last_bin_out=True).df
    mom.features("ft1").df[["Chromosome", "Start", "End"]].__eq__(ft1)


@pytest.mark.order(99999999)
def test_Momics_remove(momics_path: str):
    mom = momics.Momics(momics_path)
    mom.remove()
    assert not Path(mom.path).exists()

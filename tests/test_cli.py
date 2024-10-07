import collections
from click.testing import CliRunner
import numpy as np
import pytest
import os
import pyBigWig
import shutil
import momics.cli as cli
from pybedtools import BedTool
from momics import utils
from momics import momics
from momics import multirangequery


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture(scope="session")
def path():
    tmp_dir = os.path.join(os.getcwd(), "testCLI.mom")
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    yield tmp_dir
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)


def test_create(runner, path):
    result = runner.invoke(cli.create.create, [path])
    assert result.exit_code == 0
    assert os.path.exists(path)
    result = runner.invoke(cli.tree.tree, [path])
    assert os.path.exists(path)


def test_add_chroms(runner, path, bw3):
    chroms = utils.get_chr_lengths(bw3)
    with open("chrom_lengths.txt", "w") as f:
        for chrom, length in chroms.items():
            f.write(f"{chrom}\t{length}\n")
    result = runner.invoke(cli.add.add, ["chroms", "--file", "chrom_lengths.txt", "--genome", "S288c", path])
    assert result.exit_code == 0
    result = runner.invoke(cli.ls.ls, ["--table", "chroms", path])
    assert result.exit_code == 0
    os.remove("chrom_lengths.txt")


def test_add_tracks(runner, path, bw3):
    result = runner.invoke(cli.add.add, ["tracks", "--file", f"bw1={bw3}", "-f", f"bw2={bw3}", path])
    assert result.exit_code == 0
    result = runner.invoke(cli.ls.ls, ["--table", "tracks", path])
    assert result.exit_code == 0
    mom = momics.Momics(path)
    assert mom.tracks()["label"].__eq__(["bw1", "bw2"]).all()


def test_add_sequence(runner, path, fa3):
    result = runner.invoke(cli.add.add, ["seq", "--file", fa3, path])
    assert result.exit_code == 0
    mom = momics.Momics(path)
    assert mom.seq()["chrom"].__eq__(["I", "II", "III", "IV"]).all()


def test_add_features(runner, path):
    result = runner.invoke(cli.binnify.binnify, ["--width", "10", "-s", "10", "-o", "out.bins.bed", "-c", path])
    assert result.exit_code == 0
    result = runner.invoke(cli.add.add, ["features", "--file", "bed1=out.bins.bed", path])
    assert result.exit_code == 0
    mom = momics.Momics(path)
    assert mom.features()["label"].__eq__(["bed1"]).all()
    assert mom.features()["n"].__eq__([13001]).all()
    result = runner.invoke(cli.binnify.binnify, ["--width", "10", "-s", "10", "-o", "out.bins.bed", path])
    assert result.exit_code == 0
    result = runner.invoke(cli.add.add, ["features", "--file", "bed2=out.bins.bed", path])
    assert result.exit_code == 0
    mom = momics.Momics(path)
    assert mom.features()["label"].__eq__(["bed1", "bed2"]).all()
    assert mom.features()["n"].__eq__([13001, 13005]).all()
    result = runner.invoke(cli.tree.tree, [path])
    assert result.exit_code == 0
    result = runner.invoke(cli.ls.ls, ["--table", "features", path])
    assert result.exit_code == 0
    os.remove("out.bins.bed")


def test_remove_track(runner, path):
    result = runner.invoke(cli.remove.remove, ["--track", "bw1", path])
    assert result.exit_code == 0
    result = runner.invoke(cli.ls.ls, ["--table", "tracks", path])
    assert result.exit_code == 0
    mom = momics.Momics(path)
    assert mom.tracks()["label"].__eq__(["None", "bw2"]).all()


def test_cp_track(runner, path, bw3):
    result = runner.invoke(cli.cp.cp, ["--track", "bw2", "-o", "out.bw", path])
    assert result.exit_code == 0
    assert os.path.exists("out.bw")
    mom = momics.Momics(path)
    bed = BedTool("\n".join(["I 990 1010", "I 1990 2010"]), from_string=True)
    q = multirangequery.MultiRangeQuery(mom, bed).query_tracks()

    # pybigwig version
    res = {"bw2": collections.defaultdict(list)}
    bw = pyBigWig.open(bw3)
    for interval in bed:
        str_coord = f"{interval.chrom}:{interval.start}-{interval.end}"
        res["bw2"][str_coord] = np.array(bw.values(interval.chrom, interval.start - 1, interval.end), dtype=np.float32)
    bw.close()
    res["bw2"] = dict(res["bw2"])

    assert np.allclose(q.coverage["bw2"]["I:990-1010"], res["bw2"]["I:990-1010"], atol=1e-6)
    # assert np.allclose(q.coverage["bw2"]["I:1990-2010"], res["bw2"]["I:1990-2010"], atol=1e-6)


def test_delete(runner, path):
    result = runner.invoke(cli.delete.delete, ["-y", path])
    assert result.exit_code == 0
    assert not os.path.exists(path)

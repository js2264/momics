import collections
from click.testing import CliRunner
import numpy as np
import pyranges as pr
import pytest
import os
import pyBigWig
import shutil
from momics.cli import cli
from momics import utils
from momics.momics import Momics
from momics import query


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


def test_cli_help(runner):
    result = runner.invoke(cli.cli)
    assert result.exit_code == 0
    result = runner.invoke(cli.cli, ["--help"])
    assert result.exit_code == 0


def test_create(runner, path):
    result = runner.invoke(cli.create.create, [path])
    assert result.exit_code == 0
    assert os.path.exists(path)
    result = runner.invoke(cli.tree.tree, [path])
    assert os.path.exists(path)


def test_ingest_chroms(runner, path, bw3):
    chroms = utils.get_chr_lengths(bw3)
    print(chroms)
    with open("chrom_lengths.txt", "w") as f:
        for chrom, length in chroms.items():
            f.write(f"{chrom}\t{length}\n")
    result = runner.invoke(cli.ingest.ingest, ["chroms", "--file", "chrom_lengths.txt", "--genome", "S288c", path])
    assert result.exit_code == 0
    result = runner.invoke(cli.ls.ls, ["--table", "chroms", path])
    assert result.exit_code == 0
    mom = Momics(path)
    assert mom.chroms()["chrom"].__eq__(["I", "II", "III", "IV"]).all()
    os.remove("chrom_lengths.txt")
    result = runner.invoke(cli.tree.tree, [path])
    assert len(result.output.strip().split("\n")) == 8


def test_ingest_tracks(runner, path, bw3):
    result = runner.invoke(cli.ingest.ingest, ["tracks", "--file", f"bw1={bw3}", "-f", f"bw2={bw3}", path])
    assert result.exit_code == 0
    result = runner.invoke(cli.ls.ls, ["--table", "tracks", path])
    assert result.exit_code == 0
    mom = Momics(path)
    assert mom.tracks()["label"].__eq__(["bw1", "bw2"]).all()
    result = runner.invoke(cli.tree.tree, [path])
    assert len(result.output.strip().split("\n")) == 13


def test_ingest_sequence(runner, path, fa3):
    result = runner.invoke(cli.ingest.ingest, ["seq", "--file", fa3, path])
    assert result.exit_code == 0
    mom = Momics(path)
    assert mom.seq()["chrom"].__eq__(["I", "II", "III", "IV"]).all()
    result = runner.invoke(cli.tree.tree, [path])
    assert len(result.output.strip().split("\n")) == 16


def test_ingest_features(runner, path):
    result = runner.invoke(cli.binnify.binnify, ["--width", "10", "-s", "10", "-o", "out.bins.bed", "-c", path])
    assert result.exit_code == 0
    result = runner.invoke(cli.ingest.ingest, ["features", "--file", "bed1=out.bins.bed", path])
    assert result.exit_code == 0
    mom = Momics(path)
    assert mom.features()["label"].__eq__(["bed1"]).all()
    assert mom.features()["n"].__eq__([13001]).all()
    result = runner.invoke(cli.binnify.binnify, ["--width", "10", "-s", "10", "-o", "out.bins.bed", path])
    assert result.exit_code == 0
    result = runner.invoke(cli.ingest.ingest, ["features", "--file", "bed2=out.bins.bed", path])
    assert result.exit_code == 0
    mom = Momics(path)
    assert mom.features()["label"].__eq__(["bed1", "bed2"]).all()
    assert mom.features()["n"].__eq__([13001, 13005]).all()
    result = runner.invoke(cli.tree.tree, [path])
    assert len(result.output.strip().split("\n")) == 21
    assert result.exit_code == 0
    result = runner.invoke(cli.ls.ls, ["--table", "features", path])
    assert result.exit_code == 0


def test_query_sequence(runner, path):
    result = runner.invoke(cli.cli, ["query", "seq", "-c", "I:0-10", path])
    assert result.output == ">I:0-10\nATCGATCGAT\n"
    result = runner.invoke(cli.cli, ["query", "seq", "-f", "out.bins.bed", path])
    assert result.output[0:28] == ">I:0-10\nATCGATCGAT\n>I:10-20\n"
    result = runner.invoke(cli.cli, ["query", "seq", "-f", "out.bins.bed", "-o", "out.fa", path])
    assert result.exit_code == 0
    os.remove("out.fa")


def test_query_tracks(runner, path):
    result = runner.invoke(cli.cli, ["query", "tracks", "-c", "I:0-10", path])
    assert result.output[0:40] == "range\tchrom\tposition\tbw1\tbw2\nI:0-10\tI\t0\t"
    result = runner.invoke(cli.cli, ["query", "tracks", "-f", "out.bins.bed", path])
    assert result.output[0:40] == "range\tchrom\tposition\tbw1\tbw2\nI:0-10\tI\t0\t"
    result = runner.invoke(cli.cli, ["query", "tracks", "-f", "out.bins.bed", "-o", "out.tsv", path])
    assert result.exit_code == 0
    os.remove("out.tsv")
    os.remove("out.bins.bed")


def test_remove_track(runner, path):
    result = runner.invoke(cli.remove.remove, ["--track", "bw1", path])
    assert result.exit_code == 0
    result = runner.invoke(cli.ls.ls, ["--table", "tracks", path])
    assert result.exit_code == 0
    mom = Momics(path)
    assert mom.tracks()["label"].__eq__(["None", "bw2"]).all()


def test_cp_track(runner, path, bw3):
    result = runner.invoke(cli.cp.cp, ["--type", "track", "--label", "bw2", "-o", "out.bw", "-f", path])
    print(result.output)
    assert result.exit_code == 0
    assert os.path.exists("out.bw")
    mom = Momics(path)
    bed = pr.from_dict({"Chromosome": ["I", "I"], "Start": [990, 1990], "End": [1010, 2010]})
    q = query.MomicsQuery(mom, bed).query_tracks()

    # pybigwig version
    res = {"bw2": collections.defaultdict(list)}
    bw = pyBigWig.open(bw3)
    for _, interval in bed.df.iterrows():
        str_coord = f"{interval.Chromosome}:{interval.Start}-{interval.End}"
        res["bw2"][str_coord] = np.array(bw.values(interval.Chromosome, interval.Start, interval.End), dtype=np.float32)
    bw.close()
    res["bw2"] = dict(res["bw2"])

    assert np.allclose(q.coverage["bw2"]["I:990-1010"], res["bw2"]["I:990-1010"], atol=1e-6)
    # assert np.allclose(q.coverage["bw2"]["I:1990-2010"], res["bw2"]["I:1990-2010"], atol=1e-6)
    os.remove("out.bw")


def test_cp_features(runner, path):
    result = runner.invoke(cli.cp.cp, ["--type", "features", "--label", "bed2", "-o", "out.bed", "-f", path])
    assert result.exit_code == 0
    assert os.path.exists("out.bed")
    bed = pr.read_bed("out.bed")
    assert ("I", 300, 310) == (bed.df.iloc[30].Chromosome, bed.df.iloc[30].Start, bed.df.iloc[30].End)


def test_cp_seq(runner, path):
    result = runner.invoke(cli.cp.cp, ["--type", "sequence", "-o", "out.fa", "-f", path])
    assert result.exit_code == 0
    assert os.path.exists("out.fa")


@pytest.mark.order(3)
def test_config(runner):
    assert os.getenv("AWS_ACCESS_KEY_ID") is not None
    assert os.getenv("AWS_SECRET_ACCESS_KEY") is not None
    result = runner.invoke(cli.cli, ["config", "s3", "list"])
    assert result.exit_code == 0
    result = runner.invoke(cli.cli, ["config", "s3", "set", "aws_access_key_id", "ABCD"])
    assert result.exit_code == 0
    assert "Set aws_access_key_id to ABCD" in result.output
    result = runner.invoke(cli.cli, ["config", "s3", "get", "aws_access_key_id"])
    assert result.exit_code == 0
    assert "aws_access_key_id: ABCD" in result.output


def test_consolidate(runner, path):
    m = Momics(path)
    s = m.size()
    result = runner.invoke(cli.consolidate.consolidate, ["--vacuum", path])
    assert result.exit_code == 0
    sf = m.size()
    assert sf < s


def test_manifest(runner, path):
    result = runner.invoke(cli.manifest.manifest, ["--output", "out.json", "-f", path])
    assert result.exit_code == 0
    assert os.path.exists("out.json")
    os.remove("out.json")


def test_delete(runner, path):
    result = runner.invoke(cli.delete.delete, ["-y", "oiasudhncoaisuhmdcoiaushcd"])
    assert result.output == "Repository oiasudhncoaisuhmdcoiaushcd does not exist.\n"
    result = runner.invoke(cli.delete.delete, ["-y", path])
    assert result.exit_code == 0
    assert not os.path.exists(path)

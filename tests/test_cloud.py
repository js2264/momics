import platform
import time
import pandas as pd
import pytest
import os
import tiledb

import momics
from momics import logging
from momics import utils
from momics import config
from momics.multirangequery import MultiRangeQuery
import dotenv

dotenv.load_dotenv()


@pytest.fixture(scope="session", autouse=True)
def load_env_variables():
    dotenv.load_dotenv()


@pytest.mark.skipif(not os.getenv("GITHUB_ACTIONS"), reason="only on GHA")
@pytest.mark.order(3)
def test_s3_IO(fa1: str, bw1: str):
    assert os.getenv("AWS_ACCESS_KEY_ID") is not None
    assert os.getenv("AWS_SECRET_ACCESS_KEY") is not None

    conf = config.MomicsConfig(
        s3=config.S3Config(
            region="eu-west-3",
            access_key_id=os.getenv("AWS_ACCESS_KEY_ID"),
            secret_access_key=os.getenv("AWS_SECRET_ACCESS_KEY"),
        )
    )

    v = platform.python_version()
    s = platform.system()
    repo = "momics/pytest.momics"
    provider = "s3"
    uri = f"{provider}://{repo}-{s}-{v}"
    conf.vfs.remove_dir(uri)
    mom = momics.Momics(uri, config=conf)

    # Add chroms
    print(bw1)
    chroms = utils.get_chr_lengths(bw1)
    print(chroms)
    mom.add_chroms(chroms)
    print(mom.chroms())
    out = pd.DataFrame(
        {
            "chrom_index": [0, 1, 2],
            "chrom": ["I", "II", "III"],
            "length": [10000, 20000, 30000],
        }
    )
    assert mom.chroms().__eq__(out).all().all()

    # Add seq
    mom.add_sequence(fa1, tile=10000)

    # Add tracks
    mom.add_tracks({"bw1": bw1}, tile=10000)
    out = pd.DataFrame(
        {
            "idx": [0],
            "label": ["bw1"],
            "path": [bw1],
        }
    )
    assert mom.tracks().__eq__(out).all().all()

    # Query tracks
    q = MultiRangeQuery(mom, "I:1-10").query_tracks()
    assert len(q.coverage) == 1
    assert len(q.coverage["bw1"]["I:1-10"]) == 10
    assert q.to_df()["chrom"].__eq__(pd.Series(["I"] * 10)).all()

    # Query sequences
    q.query_sequence()
    assert len(q.seq) == 1
    assert q.seq["nucleotide"]["I:1-10"] == "ATCGATCGAT"

    ## Purge existing repo
    res = mom.remove()
    assert res is True


# @pytest.mark.skip()
@pytest.mark.skipif(not os.getenv("GITHUB_ACTIONS"), reason="only on GHA")
@pytest.mark.order(3)
def test_gcs_IO(fa1: str, bw1: str):

    conf = config.MomicsConfig(
        gcs=config.GCSConfig(
            project_id="momics-437213",
            credentials="gcs_creds.json",
        )
    )

    v = platform.python_version()
    s = platform.system()
    repo = "momics/pytest.momics"
    provider = "gcs"
    uri = f"{provider}://{repo}-{s}-{v}"
    conf.vfs.remove_dir(uri)
    mom = momics.Momics(uri, config=conf)

    # Add chroms
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

    # Add seq
    mom.add_sequence(fa1, tile=10000)

    # Add tracks
    mom.add_tracks({"bw1": bw1}, tile=10000)
    out = pd.DataFrame(
        {
            "idx": [0],
            "label": ["bw1"],
            "path": [bw1],
        }
    )
    assert mom.tracks().__eq__(out).all().all()

    # Query tracks
    q = MultiRangeQuery(mom, "I:1-10").query_tracks()
    assert len(q.coverage) == 1
    assert len(q.coverage["bw1"]["I:1-10"]) == 10
    assert q.to_df()["chrom"].__eq__(pd.Series(["I"] * 10)).all()

    # Query sequences
    q.query_sequence()
    assert len(q.seq) == 1
    assert q.seq["nucleotide"]["I:1-10"] == "ATCGATCGAT"

    ## Purge existing repo
    res = mom.remove()
    assert res is True


@pytest.mark.skipif(not os.getenv("GITHUB_ACTIONS"), reason="only on GHA")
@pytest.mark.order(3)
def test_azure_IO(fa1: str, bw1: str):
    conf = config.MomicsConfig(
        azure=config.AzureConfig(
            account_name="momics2",
            account_key=os.getenv("AZURE_STORAGE_ACCOUNT_KEY"),
        )
    )

    v = platform.python_version()
    s = platform.system()
    repo = "momics/pytest.momics"
    provider = "azure"
    uri = f"{provider}://{repo}-{s}-{v}"

    def remove_directory_until_success(vfs, dir_uri, max_retries=10, retry_delay=2):
        attempts = 0
        while attempts < max_retries:
            try:
                vfs.remove_dir(dir_uri)
                logging.logger.info(f"Successfully removed directory: {dir_uri}")
                break
            except tiledb.TileDBError as e:
                attempts += 1
                print(
                    f"Attempt {attempts}: Failed to remove directory '{dir_uri}'. Error: {e}"
                )
                if attempts < max_retries:
                    print(f"Retrying in {retry_delay} seconds...")
                    time.sleep(retry_delay)
                else:
                    print(
                        f"Max retries reached. Directory '{dir_uri}' could not be removed."
                    )
                    raise e

    remove_directory_until_success(conf.vfs, uri)

    mom = momics.Momics(uri, config=conf)

    # Add chroms
    print(bw1)
    chroms = utils.get_chr_lengths(bw1)
    print(chroms)
    mom.add_chroms(chroms)
    out = pd.DataFrame(
        {
            "chrom_index": [0, 1, 2],
            "chrom": ["I", "II", "III"],
            "length": [10000, 20000, 30000],
        }
    )
    assert mom.chroms().__eq__(out).all().all()

    # Add seq
    mom.add_sequence(fa1, threads=2, tile=10000)

    # Add tracks
    mom.add_tracks({"bw1": bw1}, threads=2, tile=10000)
    out = pd.DataFrame(
        {
            "idx": [0],
            "label": ["bw1"],
            "path": [bw1],
        }
    )
    assert mom.tracks().__eq__(out).all().all()

    # Query tracks
    q = MultiRangeQuery(mom, "I:1-10").query_tracks()
    assert len(q.coverage) == 1
    assert len(q.coverage["bw1"]["I:1-10"]) == 10
    assert q.to_df()["chrom"].__eq__(pd.Series(["I"] * 10)).all()

    # Query sequences
    q.query_sequence()
    assert len(q.seq) == 1
    assert q.seq["nucleotide"]["I:1-10"] == "ATCGATCGAT"

    ## Purge existing repo
    res = mom.remove()
    assert res is True

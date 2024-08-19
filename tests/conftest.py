import os

import numpy as np
import pyBigWig
import pytest


@pytest.fixture(scope="session")
def momics_path(tmp_path_factory):
    p = os.path.join(tmp_path_factory.getbasetemp(), "test.momics")
    return p


@pytest.fixture(scope="session")
def bw1(tmp_path_factory):
    p = os.path.join(tmp_path_factory.getbasetemp(), "bw1")
    bw = pyBigWig.open(p, "w")
    chrom_sizes = {"I": 10000, "II": 20000, "III": 30000}
    bw.addHeader(list(chrom_sizes.items()))
    for chrom, size in chrom_sizes.items():
        intervals = [(i, i + 1000, np.random.rand()) for i in range(0, size, 1000)]
        bw.addEntries(
            [chrom] * len(intervals),
            starts=[x[0] for x in intervals],
            ends=[x[1] for x in intervals],
            values=[x[2] for x in intervals],
        )
    bw.close()
    return p


@pytest.fixture(scope="session")
def bw2(tmp_path_factory):
    p = os.path.join(tmp_path_factory.getbasetemp(), "bw2")
    bw = pyBigWig.open(p, "w")
    chrom_sizes = {"I": 10000, "II": 50000, "III": 30000, "IV": 40000}
    bw.addHeader(list(chrom_sizes.items()))
    for chrom, size in chrom_sizes.items():
        intervals = [(i, i + 1000, np.random.rand()) for i in range(0, size, 1000)]
        bw.addEntries(
            [chrom] * len(intervals),
            starts=[x[0] for x in intervals],
            ends=[x[1] for x in intervals],
            values=[x[2] for x in intervals],
        )
    bw.close()
    return p

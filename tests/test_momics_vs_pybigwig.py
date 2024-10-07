import collections
import numpy as np
import pytest
import pyBigWig
from pybedtools import BedTool

import momics
from momics.multirangequery import MultiRangeQuery


@pytest.mark.order(1)
def test_match_momics_pybigwig(momics_path: str, bw1):
    mom = momics.Momics(momics_path)
    bed = BedTool("\n".join(["I 990 1010", "I 1990 2010"]), from_string=True)

    # momics version
    q = MultiRangeQuery(mom, bed).query_tracks()

    # pybigwig version
    res = {"bw2": collections.defaultdict(list)}
    bw = pyBigWig.open(bw1)
    for interval in bed:
        str_coord = f"{interval.chrom}:{interval.start}-{interval.end}"
        res["bw2"][str_coord] = np.array(bw.values(interval.chrom, interval.start - 1, interval.end), dtype=np.float32)
    bw.close()
    res["bw2"] = dict(res["bw2"])

    assert np.allclose(q.coverage["bw2"]["I:990-1010"], res["bw2"]["I:990-1010"], atol=1e-3)
    # print(q.coverage["bw2"]["I:1990-2010"])
    # print(res["bw2"]["I:1990-2010"])
    # assert np.allclose(q.coverage["bw2"]["I:1990-2010"], res["bw2"]["I:1990-2010"], atol=1e-6)

import collections
import numpy as np
import pytest
import pyBigWig

import momics
from momics.multirangequery import MultiRangeQuery


@pytest.mark.order(1)
def test_match_momics_pybigwig(momics_path: str, bw1):
    mom = momics.Momics(momics_path)
    bed = momics.utils.parse_ucsc_coordinates(["I:990-1010", "I:1990-2010"])

    # momics version
    q = MultiRangeQuery(mom, bed).query_tracks()

    # pybigwig version
    res = {"bw2": collections.defaultdict(list)}
    bw = pyBigWig.open(bw1)
    for _, interval in bed.df.iterrows():
        str_coord = f"{interval.Chromosome}:{interval.Start}-{interval.End}"
        res["bw2"][str_coord] = np.array(bw.values(interval.Chromosome, interval.Start - 1, interval.End), dtype=np.float32)
    bw.close()
    res["bw2"] = dict(res["bw2"])

    assert np.allclose(q.coverage["bw2"]["I:990-1010"], res["bw2"]["I:990-1010"], atol=1e-3)
    # print(q.coverage["bw2"]["I:1990-2010"])
    # print(res["bw2"]["I:1990-2010"])
    # assert np.allclose(q.coverage["bw2"]["I:1990-2010"], res["bw2"]["I:1990-2010"], atol=1e-6)

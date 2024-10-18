import numpy as np
import pandas as pd
import pyranges as pr
import os
import pytest

from momics import utils


@pytest.mark.order(999)
def test_utils():
    c = utils.parse_ucsc_coordinates("I:1-10")
    df = pd.DataFrame([item.split() for item in ["I 1 10"]], columns=["Chromosome", "Start", "End"]).astype(
        {"Start": int, "End": int}
    )
    gr = pr.PyRanges(df)
    assert list(c.df[0:1].Start).__eq__(list(gr.df[0:1].Start))
    assert list(c.df[0:1].End).__eq__(list(gr.df[0:1].End))
    assert list(c.df[0:1].Chromosome).__eq__(list(gr.df[0:1].Chromosome))

    c = utils.parse_ucsc_coordinates(["I:1-10", "I:2-11"])
    print(c)
    df = pd.DataFrame([item.split() for item in ["I 1 10", "I 2 11"]], columns=["Chromosome", "Start", "End"]).astype(
        {"Start": int, "End": int}
    )
    gr = pr.PyRanges(df)
    assert list(c.df[1:2].Start).__eq__(list(gr.df[1:2].Start))
    assert list(c.df[1:2].End).__eq__(list(gr.df[1:2].End))
    assert list(c.df[1:2].Chromosome).__eq__(list(gr.df[1:2].Chromosome))

    with pytest.raises(ValueError, match=r"Invalid"):
        utils.parse_ucsc_coordinates("I:1-asdc")


@pytest.mark.order(999)
def test_dict_to_bigwig():
    bw_dict = {
        "I": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "II": [11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
    }
    utils._dict_to_bigwig(bw_dict, "out.bw")
    assert os.path.exists("out.bw")
    os.remove("out.bw")


@pytest.mark.order(999)
def test_split_ranges():
    rg = pr.from_dict(
        {
            "Chromosome": ["I", "I", "I", "I", "I", "I", "I", "I", "I", "I"],
            "Start": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            "End": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
        }
    )
    train, val = utils.split_ranges(rg, 0.6, shuffle=True)
    train, val = utils.split_ranges(rg, 0.6, shuffle=False)
    assert len(train) == 6
    assert len(val) == 4
    assert train.df.iloc[0].Start == 1
    assert train.df.iloc[-1].End == 15
    assert val.df.iloc[0].Start == 7
    assert val.df.iloc[-1].End == 19


@pytest.mark.order(999)
def test_pyranges_to_bw():

    rg = pr.from_dict({"Chromosome": ["I", "I", "I"], "Start": [1, 11, 21], "End": [20, 20, 30]})
    with pytest.raises(ValueError, match="All ranges must have the same width"):
        utils.pyranges_to_bw(rg, np.array([[1], [2], [3]]), "out.bw")

    rg = pr.from_dict({"Chromosome": ["I", "I", "I"], "Start": [1, 11, 21], "End": [10, 20, 30]})

    with pytest.raises(ValueError, match=r"Length of PyRanges object must .*"):
        utils.pyranges_to_bw(rg, np.array([[1], [2]]), "out.bw")

    with pytest.raises(ValueError, match=r"All ranges must have the same width as.*"):
        utils.pyranges_to_bw(rg, np.array([[1, 2], [2, 2], [3, 2]]), "out.bw")

    utils.pyranges_to_bw(rg, np.array([[0.1] * 10, [0.2] * 10, [0.3] * 10]), "out.bw")
    assert os.path.exists("out.bw")
    os.remove("out.bw")

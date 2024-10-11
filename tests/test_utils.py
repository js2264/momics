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

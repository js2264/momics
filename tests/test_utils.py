import os
import pybedtools
import pytest

from momics import utils


@pytest.mark.order(999)
def test_utils():
    c = utils.parse_ucsc_coordinates("I:1-10")
    pbt = pybedtools.BedTool("I 1 10", from_string=True)
    assert c[0].start.__eq__(pbt[0].start)
    assert c[0].end.__eq__(pbt[0].end)
    assert c[0].chrom.__eq__(pbt[0].chrom)

    c = utils.parse_ucsc_coordinates(["I:1-10", "I:2-11"])
    print(c)
    pbt = pybedtools.BedTool("I 1 10\nI 2 11", from_string=True)
    assert c[1].start.__eq__(pbt[1].start)
    assert c[1].end.__eq__(pbt[1].end)
    assert c[1].chrom.__eq__(pbt[1].chrom)

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

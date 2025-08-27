from momics import utils as mutils
import numpy as np


def test_one_hot_encode():
    seq = "ATCGGC"
    encoded = mutils.one_hot_encode(seq)
    expected = np.array([[0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 0, 1], [0, 0, 0, 1, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])
    assert np.array_equal(encoded, expected)


def test_one_hot_decode():
    encoded = np.array([[[0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 0, 1], [0, 0, 0, 1, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]]])
    decoded = mutils.one_hot_decode(encoded)
    expected = ["ATCGGC"]
    assert decoded == expected

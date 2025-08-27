from momics import utils as mutils
from momics import nn as mnn
import numpy as np


def test_one_hot_encode():
    seq = "ATCGGC"
    encoded = mutils.one_hot_encode(seq)
    expected = np.array([[0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 0, 1], [0, 0, 0, 1, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])
    assert np.array_equal(encoded, expected)


def test_one_hot_decode():
    encoded = np.array([[[0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]]])
    decoded = mutils.one_hot_decode(encoded)
    expected = ["ATGC"]
    assert decoded == expected


def test_cor():
    y_true = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    y_pred = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    correlation = mnn.cor(y_true, y_pred)
    assert np.isclose(correlation.numpy(), 1.0)

    y_pred = np.array([[3, 2, 1], [6, 5, 4]], dtype=np.float32)
    correlation = mnn.cor(y_true, y_pred)
    assert np.isclose(correlation.numpy(), -1.0)


def test_loss_mae_cor():
    y_true = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    y_pred = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    loss = mnn.loss_mae_cor(y_true, y_pred, alpha=0)  # loss is just correlation
    assert np.isclose(loss.numpy(), 0.0, atol=1e-6)
    loss = mnn.loss_mae_cor(y_true, y_pred, alpha=1)  # loss is just mae
    assert np.isclose(loss.numpy(), 0.0, atol=1e-6)
    loss = mnn.loss_mae_cor(y_true, y_pred, alpha=0.75)  # loss is both mae (75%) and cor (25%)
    assert np.isclose(loss.numpy(), 0.0, atol=1e-6)

    y_pred = np.array([[3, 2, 1], [6, 5, 4]], dtype=np.float32)
    loss = mnn.loss_mae_cor(y_true, y_pred, alpha=0)  # loss is just correlation
    assert np.isclose(loss.numpy(), 2.0, atol=1e-6)
    loss = mnn.loss_mae_cor(y_true, y_pred, alpha=1)  # loss is just mae
    assert np.isclose(loss.numpy(), 1.333333333, atol=1e-6)
    loss = mnn.loss_mae_cor(y_true, y_pred, alpha=0.75)  # loss is both mae (75%) and cor (25%)
    assert np.isclose(loss.numpy(), 1.5, atol=1e-6)

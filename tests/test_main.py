import os

import tiledb
import numpy as np
import pytest

import momics
from momics import api

testdir = os.path.dirname(os.path.realpath(__file__))

@pytest.mark.parametrize(
    "file_path",
    [(os.path.join(testdir, "data", "test.momics"))],
)
def test_api(file_path):
    clr = momics.api.Momics(file_path)
    print(clr)

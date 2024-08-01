import os

import pytest

import momics

testdir = os.path.dirname(os.path.realpath(__file__))


@pytest.mark.parametrize(
    "fp",
    [(os.path.join(testdir, "test.momics"))],
)
def test_Momics(fp):
    mom = momics.Momics(fp)
    print(mom)

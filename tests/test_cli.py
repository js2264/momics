import os
import tempfile

import pytest
from click.testing import CliRunner

testdir = os.path.dirname(os.path.realpath(__file__))


@pytest.mark.parametrize(
    "fp",
    ["test.momics"],
)
def test_create(fp: str):
    dp = tempfile.TemporaryDirectory().name
    os.mkdir(dp)
    p = os.path.join(dp, fp)
    runner = CliRunner()
    # result = runner.invoke(create, [p])
    # assert result.exit_code == 0

import os
import pytest
import tempfile
import pandas as pd


from momics import utils


@pytest.mark.order(999)
def test_utils():
    c = utils.parse_ucsc_coordinates("I:1-10")
    assert (
        c.__eq__(pd.DataFrame({"chrom": ["I"], "start": [1], "end": [10]})).all().all()
    )

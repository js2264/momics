import os
import tempfile

import pytest

from momics import momics


@pytest.mark.order(999)
def test_export(momics_path: str):
    p = tempfile.NamedTemporaryFile().name
    mom = momics.Momics(momics_path)
    print(mom.tracks())
    mom.export_track("bw2", p)
    assert os.path.exists(p)
    os.remove(p)
    mom.export_sequence(p)
    assert os.path.exists(p)
    os.remove(p)
    mom.export_features("ft1", p)
    assert os.path.exists(p)
    os.remove(p)

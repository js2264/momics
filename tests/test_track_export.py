import os
import tempfile

import pytest

import momics
from momics.export import export_track


@pytest.mark.order(999)
def test_export(momics_path: str, bed1: str):
    p = tempfile.NamedTemporaryFile().name
    mom = momics.Momics(momics_path)
    print(mom.tracks())
    export_track(mom, "bw2", p)
    assert os.path.exists(p)
    os.remove(p)

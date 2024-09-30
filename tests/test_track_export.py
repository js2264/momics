import os
from pathlib import Path
import pytest
import tempfile

import momics
from momics import multirangequery
from momics.export import export_track


@pytest.mark.order(999)
def test_export(momics_path: str, bed1: str):
    p = tempfile.NamedTemporaryFile().name
    mom = momics.Momics(momics_path, create=False)
    print(mom.tracks())
    export_track(mom, "bw2", p)
    assert os.path.exists(p)
    os.remove(p)

import os
import tempfile

import pytest

import momics
from momics.export import export_track, export_sequence, export_features


@pytest.mark.order(999)
def test_export(momics_path: str):
    p = tempfile.NamedTemporaryFile().name
    mom = momics.Momics(momics_path)
    print(mom.tracks())
    export_track(mom, "bw2", p)
    assert os.path.exists(p)
    os.remove(p)
    export_sequence(mom, p)
    assert os.path.exists(p)
    os.remove(p)
    export_features(mom, "ft1", p)
    assert os.path.exists(p)
    os.remove(p)

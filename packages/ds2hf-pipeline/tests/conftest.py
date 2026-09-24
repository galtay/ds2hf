from __future__ import annotations

from pathlib import Path

import pytest
from ds2hf_pipeline.mskchord import source
from mini_msk_study import MEMBERS, write_tarball


@pytest.fixture
def tarball(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """The mini study as a tarball, pinned as if it were the real one."""
    path = write_tarball(tmp_path / source.TARBALL_NAME, MEMBERS)
    monkeypatch.setattr(source, "TARBALL_SHA256", source.sha256(path))
    return path

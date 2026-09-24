"""The MSK-CHORD study bundle, as cBioPortal distributes it.

cBioPortal serves each public study as one tarball. The URL carries no
version and the study is corrected in place (datahub commit eb53cc4, "Fix to
dropped variants", replaced its mutations in 2026-07), so we pin the bytes
we tested against and refuse others until a human has looked.

The same files live in the cBioPortal/datahub git repository, most of them
as Git LFS objects. The pinned tarball is byte-identical to that repository
at `DATAHUB_COMMIT`, which gives verification a second, independent source.
"""

from __future__ import annotations

import hashlib
import tarfile
from pathlib import Path

import httpx

STUDY_ID = "msk_chord_2024"
TARBALL_NAME = f"{STUDY_ID}.tar.gz"
TARBALL_URL = f"https://datahub.assets.cbioportal.org/{TARBALL_NAME}"
TARBALL_SHA256 = "d292bae945bbf3ceae22332bfd08afb8103a8c2ab77f685483d655c4298bd9d5"

DATAHUB_REPO = "cBioPortal/datahub"
DATAHUB_COMMIT = "eb53cc4a9b69fac59e3fa8db9d5204b2d25ba73e"
DATAHUB_PATH = f"public/{STUDY_ID}"

LICENSE_URL = "https://creativecommons.org/licenses/by-nc-nd/4.0/"
STUDY_URL = f"https://www.cbioportal.org/study/summary?id={STUDY_ID}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(raw_dir: Path) -> Path:
    """Download the study tarball into `raw_dir` unless the pinned bytes are there.

    Raises on any other bytes, keeping them beside the pinned name so they
    can be inspected: a new tarball means MSK or cBioPortal changed the
    study, and the pin should move only after a review.
    """
    raw_dir.mkdir(parents=True, exist_ok=True)
    out = raw_dir / TARBALL_NAME
    if out.exists() and sha256(out) == TARBALL_SHA256:
        return out
    partial = out.with_suffix(".part")
    with httpx.stream("GET", TARBALL_URL, timeout=300.0, follow_redirects=True) as response:
        response.raise_for_status()
        with partial.open("wb") as fh:
            for chunk in response.iter_bytes():
                fh.write(chunk)
    got = sha256(partial)
    if got != TARBALL_SHA256:
        unpinned = raw_dir / f"{STUDY_ID}.{got[:12]}.tar.gz"
        partial.rename(unpinned)
        raise ValueError(
            f"{TARBALL_URL} now serves sha256 {got}, pinned {TARBALL_SHA256}. "
            f"Kept as {unpinned}; review the change before moving the pin."
        )
    partial.rename(out)
    return out


def read_members(tarball: Path) -> dict[str, bytes]:
    """Every regular file in the tarball, keyed by its path inside the study dir.

    Keys drop the leading `msk_chord_2024/`, so they match datahub paths
    under `DATAHUB_PATH` (e.g. `data_mutations.txt`, `case_lists/cases_all.txt`).
    """
    prefix = f"{STUDY_ID}/"
    members = {}
    with tarfile.open(tarball, "r:gz") as tar:
        for info in tar:
            if not info.isfile():
                continue
            if not info.name.startswith(prefix):
                raise ValueError(f"{tarball}: unexpected member {info.name!r}")
            fh = tar.extractfile(info)
            assert fh is not None
            members[info.name.removeprefix(prefix)] = fh.read()
    return members

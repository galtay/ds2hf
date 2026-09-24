"""Check a built MSK-CHORD tree against the source, not against our own beliefs.

Five checks, local first:

  1. `tree_clean`         — the tree holds exactly the files a build writes,
                            one config per `data_*` file in the tarball.
  2. `source_pinned`      — the shipped tarball has the pinned sha256, and
                            LICENSE is the tarball's own LICENSE.
  3. `cells_match_source` — every config cell equals its source-text cell;
                            no row, column or header line added, dropped or
                            reordered.
  4. `source_served`      — cBioPortal serves the same tarball bytes today.
  5. `source_at_datahub`  — every tarball file is byte-identical to the
                            cBioPortal/datahub repository at the pinned
                            commit, and neither side has a file the other
                            lacks.

The source text is parsed here from scratch rather than through `convert`,
so a bug shared with the builder cannot hide itself.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from decimal import Decimal
from pathlib import Path
from typing import Any

from ds2hf_pipeline.mskchord import source

_HEADER_KEYS = (b"display_name", b"description", b"datatype", b"priority")


@dataclass
class Check:
    name: str
    passed: bool
    summary: str
    details: list[str] = field(default_factory=list)


def _configs(members: dict[str, bytes]) -> dict[str, str]:
    """{config: member} for every top-level `data_*` file."""
    return {
        m.removeprefix("data_").removesuffix(".txt").replace(".", "_"): m
        for m in members
        if m.startswith("data_") and "/" not in m
    }


def check_tree_clean(processed: Path, members: dict[str, bytes]) -> Check:
    expected = {
        Path("README.md"),
        Path("LICENSE"),
        Path("source") / source.TARBALL_NAME,
        *(Path(c) / "data.parquet" for c in _configs(members)),
    }
    present = {p.relative_to(processed) for p in processed.rglob("*") if p.is_file()}
    details = [f"  unexpected: {p}" for p in sorted(present - expected)]
    details += [f"  missing:    {p}" for p in sorted(expected - present)]
    summary = f"{len(present)} files, {len(expected)} expected"
    return Check("tree_clean", not details, summary, details)


def check_source_pinned(processed: Path, members: dict[str, bytes]) -> Check:
    got = source.sha256(processed / "source" / source.TARBALL_NAME)
    license_ok = (processed / "LICENSE").read_bytes() == members.get("LICENSE")
    details = [f"  shipped {got}", f"  pinned  {source.TARBALL_SHA256}"]
    if not license_ok:
        details.append("  LICENSE differs from the tarball's LICENSE")
    ok = got == source.TARBALL_SHA256 and license_ok
    summary = "tarball and LICENSE as pinned" if ok else "mismatch"
    return Check("source_pinned", ok, summary, details)


def _cell_matches(text: str, value: Any, kind: str) -> bool:
    if text == "":
        return value is None
    if kind == "string":
        return value == text
    if kind == "int64":
        return isinstance(value, int) and str(value) == text
    if kind == "double":
        return isinstance(value, float) and Decimal(repr(value)) == Decimal(text)
    return False


def check_cells_match_source(processed: Path, members: dict[str, bytes]) -> Check:
    """Every cell of every config equals the source text it came from.

    Text columns must equal the text; integers must print back to it;
    floats must equal it as decimal numbers, since their formatting is not
    kept. A clinical file's four `#` lines must equal each field's metadata.
    """
    import pyarrow.parquet as pq

    details, compared = [], 0
    for config, member in sorted(_configs(members).items()):
        lines = re.split(r"\r?\n", members[member].decode("utf-8"))
        if lines[-1] == "":
            lines.pop()
        n_meta = next(i for i, line in enumerate(lines) if not line.startswith("#"))
        meta = [line[1:].split("\t") for line in lines[:n_meta]]
        header, *body = (line.split("\t") for line in lines[n_meta:])

        table = pq.read_table(processed / config / "data.parquet")
        if table.column_names != header:
            details.append(f"  {config}: column names or order differ from {member}")
            continue
        if table.num_rows != len(body):
            details.append(f"  {config}: {table.num_rows:,} rows, {member} has {len(body):,}")
            continue
        for j, name in enumerate(header):
            fld = table.schema.field(name)
            if meta:
                want = {k: line[j].encode() for k, line in zip(_HEADER_KEYS, meta, strict=True)}
                if dict(fld.metadata or {}) != want:
                    details.append(f"  {config}.{name}: field metadata differs from header lines")
            kind = str(fld.type)
            values = table.column(name).to_pylist()
            bad = [
                i
                for i, (row, value) in enumerate(zip(body, values, strict=True))
                if not _cell_matches(row[j], value, kind)
            ]
            compared += len(values)
            if bad:
                i = bad[0]
                details.append(
                    f"  {config}.{name}: {len(bad):,} cell(s) differ, e.g. row {i}: "
                    f"source {body[i][j]!r} vs ours {values[i]!r}"
                )
    return Check(
        "cells_match_source",
        not details,
        f"{compared:,} cells compared across {len(_configs(members))} configs",
        details,
    )


def check_source_served(processed: Path) -> Check:
    import httpx

    response = httpx.get(source.TARBALL_URL, timeout=300.0, follow_redirects=True)
    served = hashlib.sha256(response.content).hexdigest() if response.status_code == 200 else None
    ok = served == source.TARBALL_SHA256
    return Check(
        "source_served",
        ok,
        "cBioPortal serves the pinned bytes" if ok else "cBioPortal serves different bytes",
        [f"  served {served} (HTTP {response.status_code}, {source.TARBALL_URL})"],
    )


def check_source_at_datahub(members: dict[str, bytes]) -> Check:
    """Each tarball file equals the datahub file at the pinned commit.

    Most datahub files are Git LFS pointers; their `oid` is the sha256 of
    the real bytes, so comparing hashes needs no LFS download. Small files
    are stored directly and compared byte for byte.
    """
    import httpx

    api = f"https://api.github.com/repos/{source.DATAHUB_REPO}/contents"
    raw = f"https://raw.githubusercontent.com/{source.DATAHUB_REPO}/{source.DATAHUB_COMMIT}"
    details = []
    with httpx.Client(timeout=60.0, follow_redirects=True) as client:
        listed, dirs = set(), [source.DATAHUB_PATH]
        while dirs:
            d = dirs.pop()
            resp = client.get(f"{api}/{d}", params={"ref": source.DATAHUB_COMMIT})
            resp.raise_for_status()
            for entry in resp.json():
                if entry["type"] == "dir":
                    dirs.append(entry["path"])
                else:
                    listed.add(entry["path"].removeprefix(f"{source.DATAHUB_PATH}/"))
        for path in sorted(listed ^ set(members)):
            side = "datahub" if path in listed else "tarball"
            details.append(f"  only in {side}: {path}")
        for path in sorted(listed & set(members)):
            resp = client.get(f"{raw}/{source.DATAHUB_PATH}/{path}")
            resp.raise_for_status()
            body = resp.content
            if body.startswith(b"version https://git-lfs"):
                oid = re.search(rb"oid sha256:([0-9a-f]{64})", body)
                ours = hashlib.sha256(members[path]).hexdigest()
                same = oid is not None and oid.group(1).decode() == ours
            else:
                same = body == members[path]
            if not same:
                details.append(f"  differs: {path}")
    return Check(
        "source_at_datahub",
        not details,
        f"{len(members)} files vs {source.DATAHUB_REPO}@{source.DATAHUB_COMMIT[:7]}",
        details,
    )


def verify(processed: Path, *, network: bool = True) -> list[Check]:
    members = source.read_members(processed / "source" / source.TARBALL_NAME)
    checks = [
        check_tree_clean(processed, members),
        check_source_pinned(processed, members),
        check_cells_match_source(processed, members),
    ]
    if network:
        checks += [check_source_served(processed), check_source_at_datahub(members)]
    return checks

from __future__ import annotations

import io
import tarfile
from pathlib import Path
from unittest.mock import MagicMock

import pytest
from tcga2hf_pipeline.gdc import GDCClient, and_, eq, in_


def test_filter_helpers() -> None:
    assert eq("a", "x") == {"op": "=", "content": {"field": "a", "value": "x"}}
    assert in_("a", ["x", "y"]) == {"op": "in", "content": {"field": "a", "value": ["x", "y"]}}
    combined = and_(eq("a", "x"), eq("b", "y"))
    assert combined["op"] == "and"
    assert len(combined["content"]) == 2


def _make_bulk_tar(files: dict[str, bytes]) -> bytes:
    """Build the GDC bulk-download tar.gz layout: MANIFEST.txt + <uuid>/<file_name>."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        manifest_lines = ["id\tfilename\tmd5\tsize\tstate"]
        for path, payload in files.items():
            data = io.BytesIO(payload)
            info = tarfile.TarInfo(name=path)
            info.size = len(payload)
            tar.addfile(info, data)
            uuid_part = path.split("/", 1)[0] if "/" in path else path
            manifest_lines.append(f"{uuid_part}\t{path}\tdeadbeef\t{len(payload)}\tvalidated")
        manifest_bytes = "\n".join(manifest_lines).encode()
        info = tarfile.TarInfo(name="MANIFEST.txt")
        info.size = len(manifest_bytes)
        tar.addfile(info, io.BytesIO(manifest_bytes))
    return buf.getvalue()


def test_bulk_download_extracts_to_flat_layout(tmp_path: Path) -> None:
    """bulk_download should strip the <uuid>/ prefix and write <out_dir>/<file_name>."""
    tar_bytes = _make_bulk_tar(
        {
            "uuid-a/sample_a.maf.gz": b"PAYLOAD_A",
            "uuid-b/sample_b.maf.gz": b"PAYLOAD_B",
        }
    )

    fake_resp = MagicMock()
    fake_resp.raise_for_status.return_value = None
    fake_resp.iter_bytes.return_value = iter([tar_bytes])
    fake_stream = MagicMock()
    fake_stream.__enter__.return_value = fake_resp
    fake_stream.__exit__.return_value = None

    client = GDCClient.__new__(GDCClient)
    client._client = MagicMock()
    client._client.stream.return_value = fake_stream

    client.bulk_download(["uuid-a", "uuid-b"], tmp_path)

    assert (tmp_path / "sample_a.maf.gz").read_bytes() == b"PAYLOAD_A"
    assert (tmp_path / "sample_b.maf.gz").read_bytes() == b"PAYLOAD_B"
    # MANIFEST.txt has no '/' in its name, so it should not be written.
    assert not (tmp_path / "MANIFEST.txt").exists()


def test_bulk_download_rejects_single_id(tmp_path: Path) -> None:
    """Single-id POST returns the raw file (not tar.gz); reject loudly."""
    client = GDCClient.__new__(GDCClient)
    client._client = MagicMock()
    with pytest.raises(ValueError, match="bulk_download needs >=2 ids"):
        client.bulk_download(["uuid-a"], tmp_path)


def test_bulk_download_retries_on_truncated_tar(tmp_path: Path) -> None:
    """A truncated tar.gz (GDC drops the connection mid-stream) raises EOFError
    inside tarfile; the retry decorator must redrive the batch instead of
    surfacing the error to the caller."""
    good_tar = _make_bulk_tar({"uuid-a/file_a.bin": b"AAAA", "uuid-b/file_b.bin": b"BBBB"})
    truncated = good_tar[: len(good_tar) // 2]  # cut off mid-stream

    call_count = {"n": 0}

    def fake_stream(method: str, url: str, **kwargs: object) -> MagicMock:
        call_count["n"] += 1
        body = truncated if call_count["n"] == 1 else good_tar
        resp = MagicMock()
        resp.raise_for_status.return_value = None
        resp.iter_bytes.return_value = iter([body])
        ctx = MagicMock()
        ctx.__enter__.return_value = resp
        ctx.__exit__.return_value = None
        return ctx

    client = GDCClient.__new__(GDCClient)
    client._client = MagicMock()
    client._client.stream.side_effect = fake_stream

    client.bulk_download(["uuid-a", "uuid-b"], tmp_path)

    assert call_count["n"] == 2  # retried once after the truncated response
    assert (tmp_path / "file_a.bin").read_bytes() == b"AAAA"
    assert (tmp_path / "file_b.bin").read_bytes() == b"BBBB"


def test_bulk_download_merges_trailing_singleton(tmp_path: Path) -> None:
    """A tail batch of 1 must be folded into the prior batch (no single-id POSTs)."""
    posted: list[list[str]] = []

    def fake_stream(method: str, url: str, **kwargs: object) -> MagicMock:
        ids = kwargs["json"]["ids"]
        posted.append(list(ids))
        files = {f"{u}/file_{u}.bin": f"P{u}".encode() for u in ids}
        resp = MagicMock()
        resp.raise_for_status.return_value = None
        resp.iter_bytes.return_value = iter([_make_bulk_tar(files)])
        ctx = MagicMock()
        ctx.__enter__.return_value = resp
        ctx.__exit__.return_value = None
        return ctx

    client = GDCClient.__new__(GDCClient)
    client._client = MagicMock()
    client._client.stream.side_effect = fake_stream

    # 5 ids with batch_size=2 would produce [2, 2, 1]; we expect [2, 3].
    client.bulk_download(["a", "b", "c", "d", "e"], tmp_path, batch_size=2)
    assert [len(b) for b in posted] == [2, 3]
    assert sorted(p for batch in posted for p in batch) == ["a", "b", "c", "d", "e"]
    for u in "abcde":
        assert (tmp_path / f"file_{u}.bin").read_bytes() == f"P{u}".encode()


@pytest.mark.network
def test_cases_smoke() -> None:
    """Live GDC API smoke test. Skip with `pytest -m 'not network'`."""
    with GDCClient() as c:
        cases = c.cases(
            filters=eq("project.project_id", "TCGA-CHOL"),
            fields=["case_id", "submitter_id"],
            page_size=3,
        )
    assert len(cases) >= 3
    assert all("case_id" in case for case in cases)


@pytest.mark.network
def test_bulk_download_smoke(tmp_path: Path) -> None:
    """Live GDC bulk download: pull a few open-access files and check they land."""
    from tcga2hf_pipeline.gdc import and_ as _and
    from tcga2hf_pipeline.gdc import eq as _eq

    with GDCClient() as c:
        hits = c.files(
            filters=_and(
                _eq("cases.project.project_id", "TCGA-CHOL"),
                _eq("access", "open"),
                _eq("data_type", "Masked Somatic Mutation"),
            ),
            fields=["file_id", "file_name", "file_size"],
            page_size=3,
        )
        ids = [h["file_id"] for h in hits[:3]]
        names = [h["file_name"] for h in hits[:3]]
        sizes = [h["file_size"] for h in hits[:3]]
        assert len(ids) == 3
        c.bulk_download(ids, tmp_path)

    for name, size in zip(names, sizes, strict=True):
        target = tmp_path / name
        assert target.exists(), name
        assert target.stat().st_size == size, (name, target.stat().st_size, size)


def _versions_client(records: list[dict], calls: list[list[str]]) -> GDCClient:
    """A GDCClient whose POST /files/versions returns `records`, logging id batches."""

    def post(path: str, json: dict, timeout: float) -> MagicMock:  # noqa: A002
        assert path == "/files/versions"
        calls.append(json["ids"])
        resp = MagicMock()
        resp.raise_for_status.return_value = None
        resp.json.return_value = [r for r in records if r["id"] in json["ids"]]
        return resp

    client = GDCClient.__new__(GDCClient)
    client._client = MagicMock()
    client._client.post.side_effect = post
    return client


def test_versions_maps_ids_to_records() -> None:
    records = [
        {"id": "a", "version": "1", "release": "36.0", "latest_id": "a"},
        {"id": "b", "version": "2", "release": "40.0", "latest_id": "b"},
    ]
    client = _versions_client(records, calls=[])
    out = client.versions(["a", "b"])
    assert out["a"]["release"] == "36.0"
    assert out["b"]["version"] == "2"


def test_versions_chunks_large_id_lists() -> None:
    """Version lookups must not send one unbounded POST for a big project."""
    ids = [str(i) for i in range(1250)]
    records = [{"id": i, "version": "1", "release": "36.0", "latest_id": i} for i in ids]
    calls: list[list[str]] = []
    client = _versions_client(records, calls)

    out = client.versions(ids, chunk_size=500)
    assert len(out) == 1250
    assert [len(c) for c in calls] == [500, 500, 250]


def test_versions_of_empty_list_makes_no_request() -> None:
    calls: list[list[str]] = []
    client = _versions_client([], calls)
    assert client.versions([]) == {}
    assert calls == []

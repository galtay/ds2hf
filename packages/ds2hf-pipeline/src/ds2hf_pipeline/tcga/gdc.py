from __future__ import annotations

import json
import shutil
import tarfile
import tempfile
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import httpx
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

GDC_BASE_URL = "https://api.gdc.cancer.gov"


def eq(field: str, value: str) -> dict[str, Any]:
    return {"op": "=", "content": {"field": field, "value": value}}


def in_(field: str, values: Iterable[str]) -> dict[str, Any]:
    return {"op": "in", "content": {"field": field, "value": list(values)}}


def and_(*clauses: dict[str, Any]) -> dict[str, Any]:
    return {"op": "and", "content": list(clauses)}


_RETRYABLE = (httpx.TransportError, httpx.HTTPStatusError)

# For bulk POST /data: GDC occasionally returns a truncated tar.gz at scale
# (connection drop mid-stream after a 200 OK). httpx surfaces that as
# RemoteProtocolError; gzip/tarfile catch it later as EOFError or ReadError
# during decompression. All three are transient, so retry the whole batch.
_BULK_RETRYABLE = (
    httpx.TransportError,
    httpx.HTTPStatusError,
    httpx.HTTPError,
    EOFError,
    tarfile.ReadError,
)


class GDCClient:
    def __init__(self, base_url: str = GDC_BASE_URL, timeout: float = 60.0) -> None:
        self._client = httpx.Client(base_url=base_url, timeout=timeout)

    def close(self) -> None:
        self._client.close()

    def __enter__(self) -> GDCClient:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    @retry(
        reraise=True,
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=1, max=30),
        retry=retry_if_exception_type(_RETRYABLE),
    )
    def status(self) -> dict[str, Any]:
        """Return the GDC server status dict (includes data_release, tag, commit)."""
        resp = self._client.get("/status")
        resp.raise_for_status()
        return resp.json()

    @retry(
        reraise=True,
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=1, max=60),
        retry=retry_if_exception_type(_RETRYABLE),
    )
    def dictionary(self) -> dict[str, Any]:
        """Return the GDC submission data dictionary (~10 MB; full snapshot).

        The dictionary has no explicit version field — its current contents
        are implicitly tied to whatever data release the GDC server is on
        (see `status()`). Capture both at the same moment to pin provenance.
        """
        resp = self._client.get("/v0/submission/_dictionary/_all", timeout=120.0)
        resp.raise_for_status()
        return resp.json()

    @retry(
        reraise=True,
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=1, max=30),
        retry=retry_if_exception_type(_RETRYABLE),
    )
    def _post(self, endpoint: str, payload: dict[str, Any]) -> dict[str, Any]:
        resp = self._client.post(endpoint, json=payload)
        if resp.status_code >= 500:
            resp.raise_for_status()
        resp.raise_for_status()
        return resp.json()

    def cases(
        self,
        filters: dict[str, Any],
        fields: list[str],
        expand: list[str] | None = None,
        page_size: int = 500,
    ) -> list[dict[str, Any]]:
        """Page through /cases and return all hits (raw dicts)."""
        return self._paginate("/cases", filters, fields, expand, page_size)

    def files(
        self,
        filters: dict[str, Any],
        fields: list[str],
        expand: list[str] | None = None,
        page_size: int = 500,
    ) -> list[dict[str, Any]]:
        """Page through /files and return all hits (raw dicts)."""
        return self._paginate("/files", filters, fields, expand, page_size)

    def _paginate(
        self,
        endpoint: str,
        filters: dict[str, Any],
        fields: list[str],
        expand: list[str] | None,
        page_size: int,
    ) -> list[dict[str, Any]]:
        hits: list[dict[str, Any]] = []
        offset = 0
        while True:
            payload: dict[str, Any] = {
                "filters": filters,
                "fields": ",".join(fields),
                "format": "JSON",
                "size": page_size,
                "from": offset,
            }
            if expand:
                payload["expand"] = ",".join(expand)
            data = self._post(endpoint, payload)
            page = data.get("data", {}).get("hits", [])
            hits.extend(page)
            pagination = data.get("data", {}).get("pagination", {})
            total = pagination.get("total", len(hits))
            if len(hits) >= total or not page:
                break
            offset += len(page)
        return hits

    @retry(
        reraise=True,
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=1, max=30),
        retry=retry_if_exception_type(_RETRYABLE),
    )
    @retry(
        reraise=True,
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=1, max=30),
        retry=retry_if_exception_type(_RETRYABLE),
    )
    def versions(self, file_uuids: list[str], chunk_size: int = 500) -> dict[str, dict[str, Any]]:
        """Return {file_id: version record} from POST /files/versions.

        Each record carries the file's own `version` and `release` (the GDC
        release it *first* appeared in), plus `latest_id` / `latest_version`
        for whatever supersedes it. This pins provenance per file rather
        than per fetch: the API only ever serves the current release, so the
        release we happen to fetch at says nothing about whether the bytes
        changed. A file whose `id == latest_id` is byte-identical to what
        every release since its own served, which is what actually matters
        when a modality is added to the dataset years after the rest.
        """
        out: dict[str, dict[str, Any]] = {}
        for i in range(0, len(file_uuids), chunk_size):
            chunk = file_uuids[i : i + chunk_size]
            resp = self._client.post("/files/versions", json={"ids": chunk}, timeout=120.0)
            resp.raise_for_status()
            for record in resp.json():
                if record.get("id"):
                    out[record["id"]] = record
        return out

    def download(self, file_uuid: str, out_path: Path) -> None:
        """Stream a single open-access file by UUID. Used by future genomic modules."""
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with self._client.stream("GET", f"/data/{file_uuid}") as resp:
            resp.raise_for_status()
            with out_path.open("wb") as fh:
                for chunk in resp.iter_bytes():
                    fh.write(chunk)

    def bulk_download(
        self,
        file_uuids: list[str],
        out_dir: Path,
        batch_size: int = 50,
    ) -> None:
        """Download many open-access files via batched POST /data.

        GDC returns a tar.gz of `<uuid>/<file_name>` per id (plus a MANIFEST.txt)
        when given >=2 ids; for a single id it returns the raw file. We fold any
        tail-of-1 into the prior batch so every POST has >=2 ids and we don't
        need two response handlers. Callers with exactly one file should use
        `download()` directly.

        Files are extracted as `out_dir/<file_name>` to match the layout of the
        per-file `download()` path so callers can treat the two interchangeably.
        """
        if not file_uuids:
            return
        if len(file_uuids) < 2:
            raise ValueError("bulk_download needs >=2 ids; use download() for a single file.")
        out_dir.mkdir(parents=True, exist_ok=True)

        batches = [file_uuids[i : i + batch_size] for i in range(0, len(file_uuids), batch_size)]
        # If the last batch is a single id, merge it back so every POST has >=2.
        if len(batches) >= 2 and len(batches[-1]) == 1:
            batches[-2].append(batches[-1].pop())
            batches.pop()

        for batch in batches:
            self._download_batch(batch, out_dir)

    @retry(
        reraise=True,
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=2, max=60),
        retry=retry_if_exception_type(_BULK_RETRYABLE),
    )
    def _download_batch(self, batch: list[str], out_dir: Path) -> None:
        """POST one batch to /data, stream the tar.gz response, extract files.

        Extraction happens after the full tar.gz lands in a temp file so a
        mid-stream drop surfaces as `EOFError` / `tarfile.ReadError` here
        (rather than half-corrupting an already-extracted file). Both errors
        are in `_BULK_RETRYABLE`, so the whole batch retries cleanly.
        """
        with tempfile.NamedTemporaryFile(suffix=".tar.gz", dir=out_dir) as tmp:
            with self._client.stream("POST", "/data", json={"ids": batch}, timeout=300.0) as resp:
                resp.raise_for_status()
                for chunk in resp.iter_bytes():
                    tmp.write(chunk)
            tmp.flush()
            tmp.seek(0)
            with tarfile.open(fileobj=tmp, mode="r:gz") as tar:
                for member in tar.getmembers():
                    if not member.isfile():
                        continue
                    # Tar entries are "<uuid>/<file_name>" plus a top-level
                    # MANIFEST.txt; we keep only the leaf name.
                    if "/" not in member.name:
                        continue
                    file_name = Path(member.name).name
                    src = tar.extractfile(member)
                    if src is None:
                        continue
                    target = out_dir / file_name
                    with src, target.open("wb") as dst:
                        shutil.copyfileobj(src, dst)


def write_cases_json(cases: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cases, indent=2))

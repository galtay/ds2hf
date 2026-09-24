from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from mskchord2hf_pipeline import build, card, source, verify


@pytest.fixture
def processed(tarball: Path, tmp_path: Path) -> Path:
    out = tmp_path / "processed"
    infos = build.build(tarball, out)
    card.write_card(out, infos)
    return out


def _failed(processed: Path) -> dict[str, verify.Check]:
    return {c.name: c for c in verify.verify(processed, network=False) if not c.passed}


def test_clean_build_passes(processed: Path) -> None:
    assert _failed(processed) == {}
    readme = (processed / "README.md").read_text()
    assert "license: cc-by-nc-nd-4.0" in readme
    # The wide CNA matrix ships as a file but is not a declared config.
    assert "config_name: cna\n" not in readme
    assert "config_name: cna_hg19_seg" in readme
    assert (processed / "cna" / "data.parquet").exists()


def test_build_refuses_unpinned_bytes(tarball: Path, tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(source, "TARBALL_SHA256", "0" * 64)
    with pytest.raises(ValueError, match="pinned"):
        build.build(tarball, tmp_path / "out")


def test_changed_cell_is_caught(processed: Path) -> None:
    path = processed / "timeline_labs" / "data.parquet"
    table = pq.read_table(path)
    units = table.column("LR_UNIT_MEASURE").to_pylist()
    units[0] = units[0].strip()  # trimming padding is a change, not a format
    i = table.column_names.index("LR_UNIT_MEASURE")
    pq.write_table(table.set_column(i, "LR_UNIT_MEASURE", pa.array(units)), path)

    failed = _failed(processed)
    assert set(failed) == {"cells_match_source"}
    assert "timeline_labs.LR_UNIT_MEASURE: 1 cell(s) differ" in "\n".join(
        failed["cells_match_source"].details
    )


def test_dropped_row_is_caught(processed: Path) -> None:
    path = processed / "mutations" / "data.parquet"
    pq.write_table(pq.read_table(path).slice(1), path)
    assert set(_failed(processed)) == {"cells_match_source"}


def test_changed_header_metadata_is_caught(processed: Path) -> None:
    path = processed / "clinical_patient" / "data.parquet"
    table = pq.read_table(path)
    field = table.schema.field("GENDER")
    edited = field.with_metadata({**field.metadata, b"description": b"Sex"})
    schema = table.schema.set(table.schema.get_field_index("GENDER"), edited)
    pq.write_table(table.cast(schema), path)
    assert set(_failed(processed)) == {"cells_match_source"}


def test_stray_file_is_caught(processed: Path) -> None:
    (processed / ".pytest_cache").mkdir()
    (processed / ".pytest_cache" / "README.md").write_text("x")
    assert set(_failed(processed)) == {"tree_clean"}

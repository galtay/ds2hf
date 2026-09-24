from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from ds2hf.mskchord import Chord, cna_long


def _write(root: Path, config: str, table: pa.Table) -> None:
    (root / config).mkdir(parents=True, exist_ok=True)
    pq.write_table(table, root / config / "data.parquet")


@pytest.fixture
def chord(tmp_path: Path) -> Chord:
    field = pa.field("PATIENT_ID", pa.string(), metadata={"display_name": "Patient Identifier"})
    _write(
        tmp_path,
        "clinical_patient",
        pa.table({"PATIENT_ID": ["P-1", "P-2"]}, schema=pa.schema([field])),
    )
    _write(
        tmp_path,
        "clinical_sample",
        pa.table({"SAMPLE_ID": ["S-1a", "S-1b", "S-2"], "PATIENT_ID": ["P-1", "P-1", "P-2"]}),
    )
    _write(tmp_path, "gene_panel_matrix", pa.table({"SAMPLE_ID": ["S-2", "S-1a", "S-1b"]}))
    _write(
        tmp_path,
        "mutations",
        pa.table({"Tumor_Sample_Barcode": ["S-2", "S-1a", "S-2"], "Hugo_Symbol": ["A", "B", "C"]}),
    )
    _write(tmp_path, "sv", pa.table({"Sample_Id": ["S-1b"], "Site1_Hugo_Symbol": ["ALK"]}))
    _write(tmp_path, "cna_hg19_seg", pa.table({"ID": ["S-1a"], "seg.mean": [0.1]}))
    _write(
        tmp_path,
        "cna",
        pa.table({"Hugo_Symbol": ["TP53", "MYC"], "S-1a": ["0", "NA"], "S-1b": ["-2", "-1.5"]}),
    )
    _write(
        tmp_path,
        "timeline_treatment",
        pa.table(
            {
                "PATIENT_ID": ["P-1", "P-1", "P-1", "P-1"],
                "START_DATE": pa.array([30, None, -5, 30], pa.int64()),
                "AGENT": ["x", "unknown date", "first", "y"],
            }
        ),
    )
    return Chord(tmp_path)


def test_configs_are_read_as_published(chord: Chord) -> None:
    assert "timeline_treatment" in chord.configs
    assert chord.timelines == ["timeline_treatment"]
    assert chord.table("mutations").column("Hugo_Symbol").to_pylist() == ["A", "B", "C"]
    assert chord.dictionary("clinical_patient") == [
        {"column": "PATIENT_ID", "display_name": "Patient Identifier"}
    ]


def test_patient_nests_samples_and_their_rows(chord: Chord) -> None:
    p = chord.patient("P-2")
    [sample] = p["samples"]
    assert sample["SAMPLE_ID"] == "S-2"
    # Source order is kept within a sample.
    assert [m["Hugo_Symbol"] for m in sample["mutations"]] == ["A", "C"]
    assert sample["sv"] == []
    assert p["timeline_treatment"] == []


def test_timeline_sorted_by_start_then_source_order(chord: Chord) -> None:
    agents = [e["AGENT"] for e in chord.patient("P-1")["timeline_treatment"]]
    # Ties keep source order; an event with no START_DATE goes last.
    assert agents == ["first", "x", "y", "unknown date"]


def test_nesting_keeps_every_row(chord: Chord) -> None:
    patients = list(chord.patients())
    assert [p["PATIENT_ID"] for p in patients] == ["P-1", "P-2"]
    samples = [s for p in patients for s in p["samples"]]
    assert sum(len(s["mutations"]) for s in samples) == chord.table("mutations").num_rows
    assert sum(len(p["timeline_treatment"]) for p in patients) == 4


def test_orphan_rows_are_refused(chord: Chord, tmp_path: Path) -> None:
    _write(tmp_path, "sv", pa.table({"Sample_Id": ["S-missing"], "Site1_Hugo_Symbol": ["ALK"]}))
    with pytest.raises(ValueError, match="S-missing"):
        Chord(tmp_path).patient_table()


def test_cna_long_drops_unprofiled_cells(chord: Chord) -> None:
    long = cna_long(chord.table("cna"))
    assert long.to_pylist() == [
        {"Hugo_Symbol": "TP53", "SAMPLE_ID": "S-1a", "cna": 0.0},
        {"Hugo_Symbol": "TP53", "SAMPLE_ID": "S-1b", "cna": -2.0},
        {"Hugo_Symbol": "MYC", "SAMPLE_ID": "S-1b", "cna": -1.5},
    ]
    kept = cna_long(chord.table("cna"), keep_unprofiled=True)
    assert kept.num_rows == 4
    assert kept.column("cna").null_count == 1


def test_cna_long_refuses_unknown_text() -> None:
    with pytest.raises(pa.ArrowInvalid):
        cna_long(pa.table({"Hugo_Symbol": ["TP53"], "S-1": ["amp"]}))


def test_patient_table_can_include_cna(chord: Chord) -> None:
    p = chord.patient("P-1", include_cna=True)
    assert [len(s["cna"]) for s in p["samples"]] == [1, 2]


def test_unknown_patient(chord: Chord) -> None:
    with pytest.raises(KeyError):
        chord.patient("P-9")

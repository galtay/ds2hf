from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
from tcga2hf.schema import TABULAR_TABLES
from tcga2hf_pipeline import pathology

# Smallest thing that is structurally a PDF. The pipeline never parses these
# bytes — it only carries them — so a real scanned report adds nothing to
# these tests beyond size.
_PDF_BYTES = b"%PDF-1.4\n1 0 obj\n<</Type/Catalog>>\nendobj\ntrailer\n%%EOF\n"

_CASE_UUID = "6e534dd3-81ac-4575-ab1e-02c21e10916d"
_SAMPLE_UUID = "ca1d5a74-0687-4d12-abcc-7f21c2e9fbb6"
_REPORT_UUID = "D2B18607-E16D-4570-9E96-5A7CBAFD79FC"
_FILE_NAME = f"TCGA-W5-AA2X.{_REPORT_UUID}.PDF"


def _build_synthetic_pathology_project(
    tmp_path: Path,
    *,
    with_associated_entity: bool = True,
) -> Path:
    """Synthesize a GDC-style pathology_reports dir: one PDF + manifest."""
    project_dir = tmp_path / "TCGA-XYZ"
    reports_dir = project_dir / "pathology_reports"
    reports_dir.mkdir(parents=True)
    (reports_dir / _FILE_NAME).write_bytes(_PDF_BYTES)

    entry = {
        "file_id": "synthetic-report-file-id",
        "file_name": _FILE_NAME,
        "file_size": len(_PDF_BYTES),
        "md5sum": "not-checked-at-load-time",
        "data_type": "Pathology Report",
        "cases": [{"case_id": _CASE_UUID, "submitter_id": "TCGA-W5-AA2X"}],
        "associated_entities": (
            [
                {
                    "entity_id": _SAMPLE_UUID,
                    "entity_type": "sample",
                    "entity_submitter_id": "TCGA-W5-AA2X-01A",
                }
            ]
            if with_associated_entity
            else []
        ),
        "_status": "downloaded",
    }
    (reports_dir / "manifest.json").write_text(json.dumps([entry]))
    return project_dir


def _patient_row() -> dict:
    return {
        "case_id": _CASE_UUID,
        "case_submitter_id": "TCGA-W5-AA2X",
        "samples": [
            {
                "sample_id": _SAMPLE_UUID,
                "submitter_id": "TCGA-W5-AA2X-01A",
                "pathology_report_uuid": _REPORT_UUID,
            }
        ],
    }


def test_load_for_project_carries_pdf_bytes_verbatim(tmp_path: Path) -> None:
    project_dir = _build_synthetic_pathology_project(tmp_path)
    by_case = pathology.load_for_project(project_dir)

    assert list(by_case) == [_CASE_UUID]
    (record,) = by_case[_CASE_UUID]
    # The whole point of this modality: bytes in == bytes out, no transform.
    assert record["pdf_bytes"] == _PDF_BYTES
    assert record["file_size"] == len(_PDF_BYTES)
    assert record["sample_id"] == _SAMPLE_UUID
    assert record["sample_submitter_id"] == "TCGA-W5-AA2X-01A"
    assert record["pathology_report_uuid"] == _REPORT_UUID


def test_load_for_project_missing_dir_returns_empty(tmp_path: Path) -> None:
    """No pathology_reports/ means build proceeds without reports."""
    assert pathology.load_for_project(tmp_path / "TCGA-NONE") == {}


def test_load_for_project_skips_manifest_only_entries(tmp_path: Path) -> None:
    """A manifest entry whose PDF wasn't downloaded is skipped, not faked."""
    project_dir = _build_synthetic_pathology_project(tmp_path)
    (project_dir / "pathology_reports" / _FILE_NAME).unlink()
    assert pathology.load_for_project(project_dir) == {}


def test_attach_populates_column_and_resolves_sample(tmp_path: Path) -> None:
    project_dir = _build_synthetic_pathology_project(tmp_path)
    rows = [_patient_row()]
    pathology.attach(rows, pathology.load_for_project(project_dir))

    (report,) = rows[0]["samples_pathology_report"]
    assert report["sample_id"] == _SAMPLE_UUID
    assert report["pdf_bytes"] == _PDF_BYTES


def test_attach_falls_back_to_filename_uuid(tmp_path: Path) -> None:
    """When GDC names no sample entity, the filename UUID resolves the FK.

    `sample.pathology_report_uuid` has always been in this dataset's schema,
    so a report can find its sample even if `associated_entities` is empty.
    """
    project_dir = _build_synthetic_pathology_project(tmp_path, with_associated_entity=False)
    by_case = pathology.load_for_project(project_dir)
    assert by_case[_CASE_UUID][0]["sample_id"] is None  # nothing to resolve from yet

    rows = [_patient_row()]
    pathology.attach(rows, by_case)
    (report,) = rows[0]["samples_pathology_report"]
    assert report["sample_id"] == _SAMPLE_UUID
    assert report["sample_submitter_id"] == "TCGA-W5-AA2X-01A"


def test_attach_leaves_empty_list_for_cases_without_reports(tmp_path: Path) -> None:
    project_dir = _build_synthetic_pathology_project(tmp_path)
    rows = [_patient_row(), {"case_id": "other-case", "case_submitter_id": "TCGA-XX-2"}]
    pathology.attach(rows, pathology.load_for_project(project_dir))
    assert rows[1]["samples_pathology_report"] == []


def test_report_uuid_parser_rejects_non_uuid_names() -> None:
    """A GDC rename must yield a null FK, never a wrong one."""
    assert pathology._pathology_report_uuid(_FILE_NAME) == _REPORT_UUID
    assert pathology._pathology_report_uuid("TCGA-W5-AA2X.not-a-uuid.PDF") is None
    assert pathology._pathology_report_uuid("no-dots-at-all") is None
    assert pathology._pathology_report_uuid("TCGA-W5-AA2X.PDF") is None
    # Right hyphen groups, non-hex content.
    assert (
        pathology._pathology_report_uuid("TCGA-W5-AA2X.ZZZZZZZZ-ZZZZ-ZZZZ-ZZZZ-ZZZZZZZZZZZZ.PDF")
        is None
    )


def test_records_satisfy_the_tabular_schema(tmp_path: Path) -> None:
    """Records must land in the declared parquet schema without coercion.

    `pdf_bytes` is the only binary column in the project, so this pins that
    pa.binary() accepts raw `bytes` straight from the loader.
    """
    project_dir = _build_synthetic_pathology_project(tmp_path)
    (record,) = pathology.load_for_project(project_dir)[_CASE_UUID]
    row = {**record, "case_id": _CASE_UUID, "case_submitter_id": "TCGA-W5-AA2X"}

    table = pa.Table.from_pylist([row], schema=TABULAR_TABLES["pathology_report"])
    assert table.num_rows == 1
    assert table.column("pdf_bytes").to_pylist() == [_PDF_BYTES]


def test_build_tables_only_filter_skips_unrequested_emitters(tmp_path: Path, monkeypatch) -> None:
    """`only` must not merely filter the output — it must skip the work.

    The expensive emitters re-read every MAF / STAR TSV in a project, so an
    append of one new table has to leave them uncalled or it saves nothing.
    """
    from tcga2hf_pipeline import tabular

    project_dir = _build_synthetic_pathology_project(tmp_path)
    called: list[str] = []

    def _spy(name, result):
        def inner(*args, **kwargs):
            called.append(name)
            return result

        return inner

    monkeypatch.setattr(tabular, "_expression_batches", _spy("expression", []))
    monkeypatch.setattr(tabular, "_mutations_rows", _spy("mutations", []))

    cases = [{"case_id": _CASE_UUID, "submitter_id": "TCGA-W5-AA2X", "samples": []}]
    tables = tabular.build_tables(cases, project_dir, only={"pathology_report"})

    assert set(tables) == {"pathology_report"}
    assert len(tables["pathology_report"]) == 1
    assert called == []


def test_build_tables_pulls_in_cases_for_survival_derived(tmp_path: Path) -> None:
    """survival_derived is projected off cases rows, so asking for it alone
    still has to materialize them."""
    from tcga2hf_pipeline import tabular

    project_dir = _build_synthetic_pathology_project(tmp_path)
    cases = [{"case_id": _CASE_UUID, "submitter_id": "TCGA-W5-AA2X", "samples": []}]
    tables = tabular.build_tables(cases, project_dir, only={"survival_derived"})
    assert "cases" in tables


def test_build_tables_without_filter_builds_everything(tmp_path: Path) -> None:
    from tcga2hf_pipeline import tabular

    project_dir = _build_synthetic_pathology_project(tmp_path)
    cases = [{"case_id": _CASE_UUID, "submitter_id": "TCGA-W5-AA2X", "samples": []}]
    tables = tabular.build_tables(cases, project_dir)
    for name in ("cases", "masked_somatic_mutation", "gene_expression_quantification",
                 "files", "survival_derived", "pathology_report"):
        assert name in tables

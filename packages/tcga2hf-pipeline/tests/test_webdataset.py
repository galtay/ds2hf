"""Tests for the per-patient WebDataset shard builder.

The interesting logic isn't packing bytes into a tar — it's the naming and
provenance rules, which are the whole point of this layout:

  - member names carry GDC's `data_type` and `file_id`, never our raw-dir
    shorthand or a positional index of our own;
  - samples are written contiguously and keys are dot-free, or WebDataset
    silently merges or splits patients;
  - the two checksums in `files.jsonl` mean different things, and the
    already-compressed formats must not be gzipped a second time.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import tarfile
from pathlib import Path

from tcga2hf_pipeline import webdataset as wds_mod

_CASE_UUID = "7e61e3b8-c617-4a5f-afdf-689521c9a670"
_SUBMITTER = "TCGA-3X-AAV9"
_EXPR_FILE = "a1b2c3d4-0000-4000-8000-000000000001"
_MAF_FILE = "a1b2c3d4-0000-4000-8000-000000000002"
_PDF_FILE = "a1b2c3d4-0000-4000-8000-000000000003"
_SEG_FILE_A = "a1b2c3d4-0000-4000-8000-000000000004"
_SEG_FILE_B = "a1b2c3d4-0000-4000-8000-000000000005"

_EXPR_TSV = b"gene_id\tgene_name\ttpm_unstranded\nENSG00000000003.15\tTSPAN6\t23.0141\n"
_SEG_TXT = b"GDC_Aliquot\tChromosome\tStart\tEnd\nx\t1\t3301765\t247650984\n"
_PDF_BYTES = b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\ntrailer\n%%EOF\n"


def _case() -> dict:
    return {
        "case_id": _CASE_UUID,
        "submitter_id": _SUBMITTER,
        "primary_site": "Bile duct",
        "disease_type": "Adenomas and Adenocarcinomas",
    }


def _manifest_entry(
    file_id: str, file_name: str, data: bytes, cases: list[dict] | None = None, **extra
) -> dict:
    owners = cases if cases is not None else [_case()]
    return {
        "file_id": file_id,
        "file_name": file_name,
        "file_size": len(data),
        "md5sum": hashlib.md5(data).hexdigest(),
        "access": "open",
        "gdc_version": "1",
        "gdc_first_release": "32.0",
        "gdc_superseded": False,
        "cases": [{"case_id": c["case_id"], "submitter_id": c["submitter_id"]} for c in owners],
        **extra,
    }


def _project(tmp_path: Path, cases: list[dict] | None = None) -> Path:
    """A minimal raw/<project>/ tree with five files across four data types.

    Every file is attached to every case in `cases`, so a multi-patient tree
    gives each patient the same modality mix.
    """
    root = tmp_path / "TCGA-CHOL"
    files = {
        "expression": [
            ("9a5e.rna_seq.augmented_star_gene_counts.tsv", _EXPR_FILE, _EXPR_TSV),
        ],
        # GDC serves MAFs already gzipped.
        "mutations": [
            ("00c2.wxs.aliquot_ensemble_masked.maf.gz", _MAF_FILE, gzip.compress(b"maf", 6)),
        ],
        "pathology_reports": [
            (f"{_SUBMITTER}.3A844132-F813.PDF", _PDF_FILE, _PDF_BYTES),
        ],
        # Two files of one data_type for one patient — the multiplicity case.
        "copy_number_masked": [
            ("BASIC_p_A.nocnv_grch38.seg.v2.txt", _SEG_FILE_A, _SEG_TXT),
            ("BASIC_p_B.nocnv_grch38.seg.v2.txt", _SEG_FILE_B, _SEG_TXT + b"y\t2\t1\t2\n"),
        ],
    }
    for raw_dir, entries in files.items():
        d = root / raw_dir
        d.mkdir(parents=True)
        manifest = []
        for file_name, file_id, data in entries:
            (d / file_name).write_bytes(data)
            manifest.append(_manifest_entry(file_id, file_name, data, cases=cases))
        (d / "manifest.json").write_text(json.dumps(manifest))
    return root


def _members(shard: Path) -> list[str]:
    with tarfile.open(shard) as tar:
        return [m.name for m in tar]


def _read(shard: Path, name: str) -> bytes:
    with tarfile.open(shard) as tar:
        return tar.extractfile(name).read()


def test_gdc_suffix_keeps_extension_not_filename_structure() -> None:
    """Interior components are GDC filename structure, not extension."""
    assert wds_mod._gdc_suffix("9a5e.rna_seq.augmented_star_gene_counts.tsv") == "tsv"
    assert wds_mod._gdc_suffix("BASIC_p_A.nocnv_grch38.seg.v2.txt") == "txt"
    # GDC's own compression suffix is part of the extension.
    assert wds_mod._gdc_suffix("00c2.wxs.aliquot_ensemble_masked.maf.gz") == "maf.gz"
    # GDC ships pathology reports with an uppercase suffix.
    assert wds_mod._gdc_suffix("TCGA-3X-AAV9.3A844132-F813.PDF") == "pdf"


def test_members_use_gdc_data_type_and_file_id(tmp_path: Path) -> None:
    """Member stems are snake_cased GDC data_type; multiplicity uses file_id."""
    out = tmp_path / "out"
    wds_mod.build_project([_case()], _project(tmp_path), out, "TCGA-CHOL")
    names = _members(out / "tcga-chol-000000.tar")

    assert f"{_SUBMITTER}.gene_expression_quantification.{_EXPR_FILE}.tsv.gz" in names
    assert f"{_SUBMITTER}.masked_somatic_mutation.{_MAF_FILE}.maf.gz" in names
    assert f"{_SUBMITTER}.pathology_report.{_PDF_FILE}.pdf" in names
    # Two copy-number files for one patient, disambiguated by GDC file_id and
    # never by a positional index of our own.
    for file_id in (_SEG_FILE_A, _SEG_FILE_B):
        assert f"{_SUBMITTER}.masked_copy_number_segment.{file_id}.txt.gz" in names
    assert not any(".0." in n or ".1." in n for n in names)
    # Our raw-directory shorthand must never reach a shard.
    shorthands = ("expression", "mutations", "pathology_reports", "copy_number_masked")
    assert not any(f".{s}." in n for n in names for s in shorthands)


def test_sample_key_is_dot_free_and_members_contiguous(tmp_path: Path) -> None:
    """WebDataset splits on the first dot and groups *successive* members."""
    out = tmp_path / "out"
    second = dict(_case(), case_id="0" * 36, submitter_id="TCGA-3X-AAVA")
    both = [_case(), second]
    wds_mod.build_project(both, _project(tmp_path, both), out, "TCGA-CHOL")
    names = _members(out / "tcga-chol-000000.tar")

    keys = [n.split(".", 1)[0] for n in names]
    assert all("." not in k for k in keys)
    # Each key occupies one unbroken run.
    runs = [k for i, k in enumerate(keys) if i == 0 or keys[i - 1] != k]
    assert len(runs) == len(set(runs))


def test_already_compressed_formats_are_not_regzipped(tmp_path: Path) -> None:
    out = tmp_path / "out"
    wds_mod.build_project([_case()], _project(tmp_path), out, "TCGA-CHOL")
    shard = out / "tcga-chol-000000.tar"
    names = _members(shard)

    assert not any(n.endswith(".pdf.gz") for n in names)
    assert not any(n.endswith(".maf.gz.gz") for n in names)
    # The PDF member is the GDC bytes untouched.
    assert _read(shard, f"{_SUBMITTER}.pathology_report.{_PDF_FILE}.pdf") == _PDF_BYTES


def test_files_jsonl_checksums_and_gzip_flag(tmp_path: Path) -> None:
    """`md5sum` is GDC's; `md5sum_member` is over the stored bytes."""
    out = tmp_path / "out"
    wds_mod.build_project([_case()], _project(tmp_path), out, "TCGA-CHOL")
    shard = out / "tcga-chol-000000.tar"
    records = {
        json.loads(line)["member"]: json.loads(line)
        for line in _read(shard, f"{_SUBMITTER}.files.jsonl").decode().splitlines()
    }

    expr = records[f"{_SUBMITTER}.gene_expression_quantification.{_EXPR_FILE}.tsv.gz"]
    assert expr["gzipped_by_pipeline"] is True
    assert expr["md5sum"] == hashlib.md5(_EXPR_TSV).hexdigest()
    stored = _read(shard, expr["member"])
    assert expr["md5sum_member"] == hashlib.md5(stored).hexdigest()
    assert gzip.decompress(stored) == _EXPR_TSV

    # GDC's own gzip: we store its bytes as-is, so both checksums agree.
    maf = records[f"{_SUBMITTER}.masked_somatic_mutation.{_MAF_FILE}.maf.gz"]
    assert maf["gzipped_by_pipeline"] is False
    assert maf["md5sum"] == maf["md5sum_member"]


def test_no_gzip_stores_verbatim_gdc_bytes(tmp_path: Path) -> None:
    out = tmp_path / "out"
    wds_mod.build_project([_case()], _project(tmp_path), out, "TCGA-CHOL", gzip_members=False)
    shard = out / "tcga-chol-000000.tar"
    member = f"{_SUBMITTER}.gene_expression_quantification.{_EXPR_FILE}.tsv"
    assert _read(shard, member) == _EXPR_TSV


def test_patient_with_no_open_access_files_still_gets_a_sample(tmp_path: Path) -> None:
    """Clinical-only patients keep a sample so index and shards stay aligned."""
    out = tmp_path / "out"
    orphan = dict(_case(), case_id="1" * 36, submitter_id="TCGA-3X-AAVZ")
    rows, _ = wds_mod.build_project([_case(), orphan], _project(tmp_path), out, "TCGA-CHOL")
    names = _members(out / "tcga-chol-000000.tar")

    assert "TCGA-3X-AAVZ.case.json" in names
    assert _read(out / "tcga-chol-000000.tar", "TCGA-3X-AAVZ.files.jsonl") == b""
    assert [r["n_files"] for r in rows if r["case_submitter_id"] == "TCGA-3X-AAVZ"] == [0]


def test_case_rows_stay_narrow(tmp_path: Path) -> None:
    """The `cases` table is identifiers + shard + size, not a modality matrix.

    Which data types a patient has belongs in that patient's `files.jsonl`;
    denormalising a column per data_type into it is what this pins against.
    """
    out = tmp_path / "out"
    rows, _ = wds_mod.build_project([_case()], _project(tmp_path), out, "TCGA-CHOL")
    (row,) = rows

    assert row["gdc_portal_url"] == f"https://portal.gdc.cancer.gov/cases/{_CASE_UUID}"
    assert row["shard"] == "data/TCGA-CHOL/tcga-chol-000000.tar"
    # n_files / n_bytes count GDC files; supplement slices are derived.
    assert row["n_files"] == 5
    assert row["n_bytes"] == sum(
        len(b) for b in (_EXPR_TSV, _SEG_TXT, _SEG_TXT + b"y\t2\t1\t2\n", _PDF_BYTES)
    ) + len(gzip.compress(b"maf", 6))
    assert set(row) == set(wds_mod.cases_schema().names)
    assert [c for c in row if c.startswith("n_")] == ["n_files", "n_bytes"]
    # Clinical attributes belong in the sample's case.json, not the index.
    assert "primary_site" not in row and "disease_type" not in row


def test_shards_roll_over_without_splitting_a_patient(tmp_path: Path) -> None:
    out = tmp_path / "out"
    cases = [
        dict(_case(), case_id=f"{i:036d}", submitter_id=f"TCGA-3X-AA{i:02d}") for i in range(4)
    ]
    rows, _ = wds_mod.build_project(
        cases, _project(tmp_path, cases), out, "TCGA-CHOL", shard_target_bytes=1
    )

    assert len(sorted(out.glob("*.tar"))) == len(cases)
    for row in rows:
        shard = out / Path(row["shard"]).name
        keys = {n.split(".", 1)[0] for n in _members(shard)}
        assert keys == {row["case_submitter_id"]}


def test_file_rows_are_a_table_of_contents_for_the_shard(tmp_path: Path) -> None:
    """`files.parquet` rows must match the sample's own `files.jsonl` exactly."""
    out = tmp_path / "out"
    _, members = wds_mod.build_project([_case()], _project(tmp_path), out, "TCGA-CHOL")
    shard = out / "tcga-chol-000000.tar"

    in_shard = [
        json.loads(line) for line in _read(shard, f"{_SUBMITTER}.files.jsonl").decode().splitlines()
    ]
    assert len(members) == len(in_shard)
    assert {m["member"] for m in members} == {r["member"] for r in in_shard}

    # Every member row is joinable back to its patient and shard.
    for row in members:
        assert row["case_submitter_id"] == _SUBMITTER
        assert row["shard"] == "data/TCGA-CHOL/tcga-chol-000000.tar"
        assert row["subset_of_gdc_file"] is False  # no supplements in this fixture

    # Projecting onto the schema drops nothing the schema declares.
    schema_names = set(wds_mod.files_schema().names)
    assert schema_names.issubset(set(members[0]) | {"n_rows"})


def test_file_rows_link_to_the_file_portal_page(tmp_path: Path) -> None:
    """A `files` row is about a file, so its portal link is the file's."""
    out = tmp_path / "out"
    cases, members = wds_mod.build_project([_case()], _project(tmp_path), out, "TCGA-CHOL")

    (case_row,) = cases
    assert case_row["gdc_portal_url"] == f"https://portal.gdc.cancer.gov/cases/{_CASE_UUID}"
    for row in members:
        assert row["gdc_portal_url"] == (f"https://portal.gdc.cancer.gov/files/{row['file_id']}")
    # The case is still reachable from a file row.
    assert {r["case_id"] for r in members} == {_CASE_UUID}


def test_metadata_member_names_match_their_cardinality(tmp_path: Path) -> None:
    """`case.json` is one object; `files.jsonl` is one object per line."""
    out = tmp_path / "out"
    wds_mod.build_project([_case()], _project(tmp_path), out, "TCGA-CHOL")
    shard = out / "tcga-chol-000000.tar"

    case_blob = _read(shard, f"{_SUBMITTER}.case.json")
    assert json.loads(case_blob)["submitter_id"] == _SUBMITTER  # a single object

    lines = _read(shard, f"{_SUBMITTER}.files.jsonl").decode().splitlines()
    assert len(lines) == 5  # one per member holding a file
    assert all(json.loads(line)["member"].startswith(_SUBMITTER) for line in lines)

    names = _members(shard)
    assert not any(".cases.json" in n for n in names)  # plural was wrong: one case

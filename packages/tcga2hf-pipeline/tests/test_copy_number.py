"""Tests for the two copy-number segment tables.

The interesting logic is FK resolution, not parsing: allele-specific files
name two aliquots and we have to pick the tumour one correctly, while masked
files name one. Both identify their subject through the file's own
`GDC_Aliquot` column rather than by reading sample-type digits out of a
barcode, and these tests pin that down.
"""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
from tcga2hf.schema import TABULAR_TABLES
from tcga2hf_pipeline import tabular

_CASE_UUID = "5fb3affa-3661-48e8-9d0f-8f0f7f6b0f11"
_SAMPLE_UUID = "e235e45d-9cbe-4c11-9a3a-1cf1d1e0aa10"
_NORMAL_SAMPLE_UUID = "aa1d5a74-0687-4d12-abcc-7f21c2e9fb00"
_TUMOR_ALIQUOT = "b3c68473-78a7-44d2-95bf-7997a7e6c0e8"
_NORMAL_ALIQUOT = "451f8cc6-05fe-48f6-9f2e-1b8c0d3a77aa"

_ALLELE_SPECIFIC_TSV = (
    "GDC_Aliquot\tChromosome\tStart\tEnd\tCopy_Number\tMajor_Copy_Number\tMinor_Copy_Number\n"
    f"{_TUMOR_ALIQUOT}\tchr1\t62920\t33086177\t1\t1\t0\n"
    f"{_TUMOR_ALIQUOT}\tchr1\t33086544\t64924592\t3\t2\t1\n"
)

_MASKED_TSV = (
    "GDC_Aliquot\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
    f"{_TUMOR_ALIQUOT}\t1\t3301765\t6028611\t1917\t-0.2829\n"
    f"{_TUMOR_ALIQUOT}\tX\t6040795\t6041047\t2\t-2.0991\n"
)


def _case() -> dict:
    """A GDC case dict with one tumour sample and one matched normal sample."""
    return {
        "case_id": _CASE_UUID,
        "submitter_id": "TCGA-W5-AA33",
        "samples": [
            {
                "sample_id": _SAMPLE_UUID,
                "submitter_id": "TCGA-W5-AA33-01A",
                "sample_type": "Primary Tumor",
                "portions": [
                    {
                        "portion_id": "portion-tumor",
                        "submitter_id": "TCGA-W5-AA33-01A-11",
                        "analytes": [
                            {
                                "aliquots": [
                                    {
                                        "aliquot_id": _TUMOR_ALIQUOT,
                                        "submitter_id": "TCGA-W5-AA33-01A-11D-A416-01",
                                    }
                                ]
                            }
                        ],
                    }
                ],
            },
            {
                "sample_id": _NORMAL_SAMPLE_UUID,
                "submitter_id": "TCGA-W5-AA33-10A",
                "sample_type": "Blood Derived Normal",
                "portions": [
                    {
                        "portion_id": "portion-normal",
                        "submitter_id": "TCGA-W5-AA33-10A-01",
                        "analytes": [
                            {
                                "aliquots": [
                                    {
                                        "aliquot_id": _NORMAL_ALIQUOT,
                                        "submitter_id": "TCGA-W5-AA33-10A-01D-A419-01",
                                    }
                                ]
                            }
                        ],
                    }
                ],
            },
        ],
    }


def _build_project(
    tmp_path: Path,
    *,
    modality: str,
    file_name: str,
    content: str,
    entities: list[dict],
    workflow_type: str,
    experimental_strategy: str | None = None,
    status: str = "downloaded",
    write_file: bool = True,
) -> Path:
    project_dir = tmp_path / "TCGA-XYZ"
    mod_dir = project_dir / modality
    mod_dir.mkdir(parents=True)
    if write_file:
        (mod_dir / file_name).write_text(content)
    entry = {
        "file_id": f"synthetic-{modality}-file-id",
        "file_name": file_name,
        "workflow_type": workflow_type,
        "experimental_strategy": experimental_strategy,
        "cases": [{"case_id": _CASE_UUID, "submitter_id": "TCGA-W5-AA33"}],
        "associated_entities": entities,
        "_status": status,
    }
    (mod_dir / "manifest.json").write_text(json.dumps([entry]))
    return project_dir


def _paired_entities() -> list[dict]:
    """Both aliquots, normal listed first — order must not decide the tumour."""
    return [
        {
            "entity_id": _NORMAL_ALIQUOT,
            "entity_type": "aliquot",
            "entity_submitter_id": "TCGA-W5-AA33-10A-01D-A419-01",
        },
        {
            "entity_id": _TUMOR_ALIQUOT,
            "entity_type": "aliquot",
            "entity_submitter_id": "TCGA-W5-AA33-01A-11D-A416-01",
        },
    ]


def _allele_specific_project(tmp_path: Path, **kwargs) -> Path:
    return _build_project(
        tmp_path,
        modality="copy_number_allele_specific",
        file_name="TCGA-CHOL.aliquot.ascat3.allelic_specific.seg.txt",
        content=_ALLELE_SPECIFIC_TSV,
        entities=_paired_entities(),
        workflow_type="ASCAT3",
        experimental_strategy="Genotyping Array",
        **kwargs,
    )


def _masked_project(tmp_path: Path, **kwargs) -> Path:
    return _build_project(
        tmp_path,
        modality="copy_number_masked",
        file_name="SPIKY_p_TCGAb.nocnv_grch38.seg.v2.txt",
        content=_MASKED_TSV,
        entities=[
            {
                "entity_id": _TUMOR_ALIQUOT,
                "entity_type": "aliquot",
                "entity_submitter_id": "TCGA-W5-AA33-01A-11D-A416-01",
            }
        ],
        workflow_type="DNAcopy",
        experimental_strategy="Genotyping Array",
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Allele-specific
# ---------------------------------------------------------------------------


def test_tumor_aliquot_comes_from_the_file_not_the_entity_order(tmp_path: Path) -> None:
    """The normal is listed first in `associated_entities`; the tumour still wins.

    This is the whole point of reading `GDC_Aliquot`: entity order is not a
    contract, and inferring tumour-vs-normal from barcode digits would be.
    """
    project = _allele_specific_project(tmp_path)
    rows = tabular._allele_specific_copy_number_segment_rows([_case()], project)

    assert len(rows) == 2
    for row in rows:
        assert row["aliquot_id"] == _TUMOR_ALIQUOT
        assert row["aliquot_submitter_id"] == "TCGA-W5-AA33-01A-11D-A416-01"
        assert row["matched_normal_aliquot_id"] == _NORMAL_ALIQUOT
        assert row["matched_normal_aliquot_submitter_id"] == "TCGA-W5-AA33-10A-01D-A419-01"
        # The sample FK must follow the tumour aliquot, not the normal one.
        assert row["sample_id"] == _SAMPLE_UUID
        assert row["sample_type"] == "Primary Tumor"


def test_allele_specific_carries_values_and_workflow_verbatim(tmp_path: Path) -> None:
    project = _allele_specific_project(tmp_path)
    rows = tabular._allele_specific_copy_number_segment_rows([_case()], project)

    first, second = rows
    assert (first["chromosome"], first["start"], first["end"]) == ("chr1", 62920, 33086177)
    assert (first["copy_number"], first["major_copy_number"], first["minor_copy_number"]) == (
        1,
        1,
        0,
    )
    assert second["copy_number"] == second["major_copy_number"] + second["minor_copy_number"]
    assert first["workflow_type"] == "ASCAT3"
    assert first["experimental_strategy"] == "Genotyping Array"


def test_allele_specific_keeps_the_chr_prefix(tmp_path: Path) -> None:
    """`chr1` here vs bare `1` in the masked table — both as the source writes them."""
    project = _allele_specific_project(tmp_path)
    rows = tabular._allele_specific_copy_number_segment_rows([_case()], project)
    assert {r["chromosome"] for r in rows} == {"chr1"}


def test_ambiguous_gdc_aliquot_column_is_skipped_not_guessed(tmp_path: Path) -> None:
    """A file naming two different aliquots has no defensible tumour FK."""
    mixed = _ALLELE_SPECIFIC_TSV + f"{_NORMAL_ALIQUOT}\tchr2\t1\t100\t2\t1\t1\n"
    project = _build_project(
        tmp_path,
        modality="copy_number_allele_specific",
        file_name="ambiguous.seg.txt",
        content=mixed,
        entities=_paired_entities(),
        workflow_type="ASCAT3",
    )
    assert tabular._allele_specific_copy_number_segment_rows([_case()], project) == []


def test_missing_matched_normal_leaves_null_rather_than_a_guess(tmp_path: Path) -> None:
    project = _build_project(
        tmp_path,
        modality="copy_number_allele_specific",
        file_name="unpaired.seg.txt",
        content=_ALLELE_SPECIFIC_TSV,
        entities=[
            {
                "entity_id": _TUMOR_ALIQUOT,
                "entity_type": "aliquot",
                "entity_submitter_id": "TCGA-W5-AA33-01A-11D-A416-01",
            }
        ],
        workflow_type="AscatNGS",
    )
    rows = tabular._allele_specific_copy_number_segment_rows([_case()], project)
    assert rows
    assert all(r["matched_normal_aliquot_id"] is None for r in rows)
    assert all(r["aliquot_id"] == _TUMOR_ALIQUOT for r in rows)


# ---------------------------------------------------------------------------
# Masked
# ---------------------------------------------------------------------------


def test_masked_segments_keep_bare_chromosome_names_as_strings(tmp_path: Path) -> None:
    """`1` must stay the string "1" — not the integer 1, and X must survive."""
    project = _masked_project(tmp_path)
    rows = tabular._masked_copy_number_segment_rows([_case()], project)

    assert [r["chromosome"] for r in rows] == ["1", "X"]
    assert all(isinstance(r["chromosome"], str) for r in rows)


def test_masked_carries_probe_counts_and_segment_means(tmp_path: Path) -> None:
    project = _masked_project(tmp_path)
    rows = tabular._masked_copy_number_segment_rows([_case()], project)

    assert (rows[0]["num_probes"], rows[0]["segment_mean"]) == (1917, -0.2829)
    assert rows[0]["aliquot_id"] == _TUMOR_ALIQUOT
    assert rows[0]["workflow_type"] == "DNAcopy"


# ---------------------------------------------------------------------------
# Shared contracts
# ---------------------------------------------------------------------------


def test_manifest_only_entries_are_skipped(tmp_path: Path) -> None:
    """A capped fetch lists files it didn't download; they have no bytes to read."""
    project = _allele_specific_project(tmp_path, status="manifest_only", write_file=False)
    assert tabular._allele_specific_copy_number_segment_rows([_case()], project) == []


def test_missing_modality_dir_yields_no_rows(tmp_path: Path) -> None:
    project = tmp_path / "TCGA-EMPTY"
    project.mkdir()
    assert tabular._allele_specific_copy_number_segment_rows([_case()], project) == []
    assert tabular._masked_copy_number_segment_rows([_case()], project) == []


def test_unknown_case_is_skipped(tmp_path: Path) -> None:
    project = _masked_project(tmp_path)
    assert tabular._masked_copy_number_segment_rows([], project) == []


def test_rows_satisfy_the_tabular_schemas(tmp_path: Path) -> None:
    """Emitted rows must load into the published schemas without coercion."""
    allele = tabular._allele_specific_copy_number_segment_rows(
        [_case()], _allele_specific_project(tmp_path / "a")
    )
    masked = tabular._masked_copy_number_segment_rows([_case()], _masked_project(tmp_path / "m"))

    for rows, table in (
        (allele, "allele_specific_copy_number_segment"),
        (masked, "masked_copy_number_segment"),
    ):
        built = pa.Table.from_pylist(rows, schema=TABULAR_TABLES[table])
        assert built.num_rows == len(rows)
        assert built.schema.equals(TABULAR_TABLES[table])


# ---------------------------------------------------------------------------
# Unmasked segments (`copy_number_segment`) — two workflows, one table
#
# DNAcopy is single-aliquot with a UUID in `GDC_Aliquot`; GATK4 CNV is paired
# with a *barcode* in `GDC_Aliquot_ID`. The barcode is the only thing that
# distinguishes tumour from matched normal in a GATK4 file, and it has to be
# matched against `entity_submitter_id` rather than parsed for sample-type
# digits — that is what these tests pin.
# ---------------------------------------------------------------------------

_UNMASKED_DNACOPY_TSV = (
    "GDC_Aliquot\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
    f"{_TUMOR_ALIQUOT}\t1\t62920\t668210\t21\t0.2874\n"
    f"{_TUMOR_ALIQUOT}\t1\t771719\t2853893\t522\t-0.1598\n"
)

_GATK4_TSV = (
    "GDC_Aliquot_ID\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
    "TCGA-W5-AA33-01A-11D-A416-01\tchr1\t17001\t828000\t59\t0.206868\n"
    "TCGA-W5-AA33-01A-11D-A416-01\tchr1\t828001\t4846000\t3809\t-0.336467\n"
)


def _dnacopy_project(tmp_path: Path, **kwargs) -> Path:
    return _build_project(
        tmp_path,
        modality="copy_number_segment_dnacopy",
        file_name="NULLS_p_TCGA.grch38.seg.v2.txt",
        content=_UNMASKED_DNACOPY_TSV,
        entities=[
            {
                "entity_id": _TUMOR_ALIQUOT,
                "entity_type": "aliquot",
                "entity_submitter_id": "TCGA-W5-AA33-01A-11D-A416-01",
            }
        ],
        workflow_type="DNAcopy",
        experimental_strategy="Genotyping Array",
        **kwargs,
    )


def _gatk4_project(tmp_path: Path, **kwargs) -> Path:
    return _build_project(
        tmp_path,
        modality="copy_number_segment_gatk4",
        file_name="uuid_wgs_gdc_realn.cr.igv.reheader.seg.txt",
        content=_GATK4_TSV,
        entities=_paired_entities(),
        workflow_type="GATK4 CNV",
        experimental_strategy="WGS",
        **kwargs,
    )


def test_unmasked_dnacopy_rows_resolve_from_gdc_aliquot(tmp_path: Path) -> None:
    rows = tabular._copy_number_segment_rows([_case()], _dnacopy_project(tmp_path))
    assert len(rows) == 2
    first = rows[0]
    assert first["aliquot_id"] == _TUMOR_ALIQUOT
    assert first["sample_id"] == _SAMPLE_UUID
    assert first["workflow_type"] == "DNAcopy"
    assert first["experimental_strategy"] == "Genotyping Array"
    # Single-aliquot workflow: no matched normal to report.
    assert first["matched_normal_aliquot_id"] is None
    assert first["matched_normal_aliquot_submitter_id"] is None
    # Bare chromosome names in this workflow, carried as written.
    assert first["chromosome"] == "1"
    assert (first["num_probes"], first["segment_mean"]) == (21, 0.2874)


def test_gatk4_tumor_resolved_by_barcode_not_entity_order(tmp_path: Path) -> None:
    rows = tabular._copy_number_segment_rows([_case()], _gatk4_project(tmp_path))
    assert len(rows) == 2
    first = rows[0]
    # The normal is listed first in associated_entities; the barcode in the
    # file's own column is what decides, so the tumour must still win.
    assert first["aliquot_id"] == _TUMOR_ALIQUOT
    assert first["matched_normal_aliquot_id"] == _NORMAL_ALIQUOT
    assert first["matched_normal_aliquot_submitter_id"] == "TCGA-W5-AA33-10A-01D-A419-01"
    assert first["workflow_type"] == "GATK4 CNV"
    assert first["experimental_strategy"] == "WGS"
    # `chr`-prefixed in this workflow, unlike DNAcopy above.
    assert first["chromosome"] == "chr1"


def test_gatk4_file_whose_barcode_matches_no_entity_is_skipped(tmp_path: Path) -> None:
    project = _build_project(
        tmp_path,
        modality="copy_number_segment_gatk4",
        file_name="uuid_wgs_gdc_realn.cr.igv.reheader.seg.txt",
        content=_GATK4_TSV.replace("TCGA-W5-AA33-01A-11D-A416-01", "TCGA-ZZ-9999-01A-11D-XXXX-01"),
        entities=_paired_entities(),
        workflow_type="GATK4 CNV",
    )
    assert tabular._copy_number_segment_rows([_case()], project) == []


def test_both_unmasked_workflows_union_into_one_table(tmp_path: Path) -> None:
    project = _dnacopy_project(tmp_path)
    _gatk4_project(tmp_path)  # same project dir, second modality
    rows = tabular._copy_number_segment_rows([_case()], project)
    assert len(rows) == 4
    assert {r["workflow_type"] for r in rows} == {"DNAcopy", "GATK4 CNV"}
    pa.Table.from_pylist(rows, schema=TABULAR_TABLES["copy_number_segment"])


# ---------------------------------------------------------------------------
# Gene-level copy number
#
# The three ASCAT workflows are paired but their TSVs carry no aliquot column,
# so the tumour is recovered from the matching allele-specific segment file.
# ABSOLUTE LiftOver names one aliquot and needs no such lookup.
# ---------------------------------------------------------------------------

_GENE_LEVEL_TSV = (
    "gene_id\tgene_name\tchromosome\tstart\tend\tcopy_number\tmin_copy_number\tmax_copy_number\n"
    "ENSG00000223972.5\tDDX11L1\tchr1\t11869\t14409\t2\t2\t2\n"
    "ENSG00000227232.5\tWASH7P\tchr1\t14404\t29570\t3\t2\t3\n"
    "ENSG00000278267.1\tMIR6859-1\tchr1\t17369\t17436\t\t\t\n"
)


def _gene_level_project(tmp_path: Path, *, workflow: str, entities: list[dict]) -> Path:
    return _build_project(
        tmp_path,
        modality="gene_level_copy_number",
        file_name=f"TCGA-CHOL.aliquot.{workflow}.gene_level_copy_number.v36.tsv",
        content=_GENE_LEVEL_TSV,
        entities=entities,
        workflow_type=workflow,
        experimental_strategy="Genotyping Array",
    )


def _gene_level_rows(cases: list[dict], project: Path) -> list[dict]:
    batches = list(tabular._gene_level_copy_number_batches(cases, project))
    return [row for b in batches for row in b.to_pylist()]


def test_absolute_gene_level_needs_no_pair_lookup(tmp_path: Path) -> None:
    """ABSOLUTE ships no segment file, so it must resolve on its own."""
    project = _gene_level_project(
        tmp_path,
        workflow="ABSOLUTE LiftOver",
        entities=[
            {
                "entity_id": _TUMOR_ALIQUOT,
                "entity_type": "aliquot",
                "entity_submitter_id": "TCGA-W5-AA33-01A-11D-A416-01",
            }
        ],
    )
    rows = _gene_level_rows([_case()], project)
    assert len(rows) == 3
    assert rows[0]["aliquot_id"] == _TUMOR_ALIQUOT
    assert rows[0]["sample_id"] == _SAMPLE_UUID
    assert rows[0]["matched_normal_aliquot_id"] is None
    assert rows[0]["gene_id"] == "ENSG00000223972.5"
    assert (rows[1]["copy_number"], rows[1]["min_copy_number"]) == (3, 2)
    # Uncalled genes stay null rather than becoming 0.
    assert rows[2]["copy_number"] is None
    # gene_name / coordinates live in `gene_model`, not here.
    assert "gene_name" not in rows[0]


def test_paired_gene_level_resolves_tumor_via_segment_file(tmp_path: Path) -> None:
    _allele_specific_project(tmp_path)  # ASCAT3 seg file naming the tumour
    project = _gene_level_project(
        tmp_path, workflow="ASCAT3", entities=_paired_entities()
    )
    rows = _gene_level_rows([_case()], project)
    assert len(rows) == 3
    assert rows[0]["aliquot_id"] == _TUMOR_ALIQUOT
    assert rows[0]["matched_normal_aliquot_id"] == _NORMAL_ALIQUOT
    assert rows[0]["workflow_type"] == "ASCAT3"


def test_paired_gene_level_is_skipped_when_the_pair_is_unresolvable(tmp_path: Path) -> None:
    """No allele-specific file to consult, so the tumour is unknown — skip.

    Guessing from entity order or barcode digits would silently mislabel
    every matched normal as the tumour.
    """
    project = _gene_level_project(
        tmp_path, workflow="ASCAT3", entities=_paired_entities()
    )
    assert _gene_level_rows([_case()], project) == []


def test_gene_level_batches_conform_to_the_published_schema(tmp_path: Path) -> None:
    _allele_specific_project(tmp_path)
    project = _gene_level_project(tmp_path, workflow="ASCAT3", entities=_paired_entities())
    batches = list(tabular._gene_level_copy_number_batches([_case()], project))
    assert batches
    for batch in batches:
        assert batch.schema.equals(TABULAR_TABLES["gene_level_copy_number"])


def test_streaming_writer_round_trips_gene_level_rows(tmp_path: Path) -> None:
    """The batch path must produce the same parquet the list path would."""
    import pyarrow.parquet as pq

    _allele_specific_project(tmp_path)
    project = _gene_level_project(tmp_path, workflow="ASCAT3", entities=_paired_entities())
    tables = {
        "gene_level_copy_number": tabular._gene_level_copy_number_batches([_case()], project)
    }
    out = tabular.write_tables(tables, tmp_path / "processed", "TCGA-XYZ")
    written = pq.read_table(out["gene_level_copy_number"])
    assert written.num_rows == 3
    assert written.schema.equals(TABULAR_TABLES["gene_level_copy_number"])
    assert written.column("gene_id").to_pylist()[0] == "ENSG00000223972.5"


def test_modality_with_no_files_writes_no_parquet(tmp_path: Path) -> None:
    tables = {"gene_level_copy_number": iter(())}
    out = tabular.write_tables(tables, tmp_path / "processed", "TCGA-XYZ")
    assert "gene_level_copy_number" not in out
    assert not (tmp_path / "processed" / "TCGA-XYZ" / "gene_level_copy_number").exists()


# ---------------------------------------------------------------------------
# Gene model — the reference table both per-gene tables join against
# ---------------------------------------------------------------------------

_EXPRESSION_TSV = (
    "# gene-model: GENCODE v36\n"
    "gene_id\tgene_name\tgene_type\tunstranded\ttpm_unstranded\n"
    "N_unmapped\t\t\t100\t0\n"
    "ENSG00000223972.5\tDDX11L1\ttranscribed_unprocessed_pseudogene\t3\t0.1\n"
    "ENSG00000227232.5\tWASH7P\tunprocessed_pseudogene\t50\t2.5\n"
    "ENSG00000198695.2\tMT-ND6\tprotein_coding\t900\t60.0\n"
)


def test_gene_model_merges_both_gdc_sources(tmp_path: Path) -> None:
    project = _build_project(
        tmp_path,
        modality="expression",
        file_name="star.rna_seq.augmented_star_gene_counts.tsv",
        content=_EXPRESSION_TSV,
        entities=[{"entity_id": _TUMOR_ALIQUOT, "entity_type": "aliquot"}],
        workflow_type="STAR - Counts",
    )
    _gene_level_project(
        tmp_path,
        workflow="ABSOLUTE LiftOver",
        entities=[{"entity_id": _TUMOR_ALIQUOT, "entity_type": "aliquot"}],
    )
    rows = tabular._gene_model_rows(project)
    by_id = {r["gene_id"]: r for r in rows}

    # The N_* alignment-summary rows are not genes.
    assert not any(g.startswith("N_") for g in by_id)
    # Union of both sources: 3 expression genes + 1 gene-level-only gene.
    assert set(by_id) == {
        "ENSG00000223972.5",
        "ENSG00000227232.5",
        "ENSG00000278267.1",
        "ENSG00000198695.2",
    }

    # A gene in both sources gets both halves.
    both = by_id["ENSG00000223972.5"]
    assert both["gene_name"] == "DDX11L1"
    assert both["gene_type"] == "transcribed_unprocessed_pseudogene"
    assert (both["chromosome"], both["start"], both["end"]) == ("chr1", 11869, 14409)

    # chrM genes are expression-only: named and typed, but no coordinates
    # invented from outside the GDC.
    mt = by_id["ENSG00000198695.2"]
    assert mt["gene_name"] == "MT-ND6"
    assert mt["gene_type"] == "protein_coding"
    assert (mt["chromosome"], mt["start"], mt["end"]) == (None, None, None)

    # A gene only the copy number files carry has coordinates but no type.
    cnv_only = by_id["ENSG00000278267.1"]
    assert cnv_only["gene_type"] is None
    assert cnv_only["chromosome"] == "chr1"

    assert [r["gene_id"] for r in rows] == sorted(by_id)
    pa.Table.from_pylist(rows, schema=TABULAR_TABLES["gene_model"])


def test_gene_model_is_empty_without_either_source(tmp_path: Path) -> None:
    project = tmp_path / "TCGA-XYZ"
    project.mkdir()
    assert tabular._gene_model_rows(project) == []

from __future__ import annotations

import pyarrow as pa
import pytest
from mini_msk_study import CLINICAL_PATIENT, CNA, MUTATIONS, SEG, TIMELINE_LABS
from mskchord2hf_pipeline import convert


def _read(member: str, text: str) -> pa.Table:
    return convert.read_table(member, text.encode())


def test_config_names() -> None:
    assert convert.config_name("data_timeline_ca_15-3_labs.txt") == "timeline_ca_15-3_labs"
    assert convert.config_name("data_cna_hg19.seg") == "cna_hg19_seg"
    assert convert.config_name("meta_cna.txt") is None
    assert convert.config_name("case_lists/cases_all.txt") is None


def test_line_endings_are_format_but_padding_is_data() -> None:
    table = _read("data_timeline_labs.txt", TIMELINE_LABS)
    # CRLF and the missing final newline end lines; neither reaches a value.
    assert table.num_rows == 3
    assert table.column("LR_UNIT_MEASURE").to_pylist() == ["ng/ml     "] * 3
    # Empty cells are null; row order is the file's.
    assert table.column("STOP_DATE").to_pylist() == [None, None, None]
    assert table.column("START_DATE").to_pylist() == [30, -5, 0]


def test_numbers_only_when_lossless() -> None:
    table = _read("data_timeline_labs.txt", TIMELINE_LABS)
    assert table.schema.field("START_DATE").type == pa.int64()
    # "2" and "0" in a column with "1.25" are floats; their text is not kept.
    assert table.schema.field("RESULT").type == pa.float64()
    assert table.column("RESULT").to_pylist() == [2.0, 1.25, 0.0]

    seg = _read("data_cna_hg19.seg", SEG)
    assert seg.column("seg.mean").to_pylist() == [-0.0004, 0.25]
    # One "X" keeps a chromosome column text; "-" keeps an allele column text.
    assert seg.schema.field("chrom").type == pa.string()
    assert _read("data_mutations.txt", MUTATIONS).schema.field("Chromosome").type == pa.string()


@pytest.mark.parametrize(
    "value",
    [
        "007",  # zero padding would be lost
        ".5",  # not the canonical form we accept
        "0.12345678901234567891",  # more digits than float64 holds
        "1_000",  # Python's float() accepts it; it is not a number in a TSV
        "nan",
    ],
)
def test_lossy_or_odd_numbers_stay_text(value: str) -> None:
    table = _read("data_x.txt", f"A\tB\n{value}\t1\n")
    assert table.schema.field("A").type == pa.string()
    assert table.column("A").to_pylist() == [value]


def test_clinical_header_lines_become_field_metadata() -> None:
    table = _read("data_clinical_patient.txt", CLINICAL_PATIENT)
    assert table.column_names == ["PATIENT_ID", "GENDER", "OS_MONTHS", "GLEASON"]
    assert table.schema.field("PATIENT_ID").metadata == {
        b"display_name": b"Patient Identifier",
        b"description": b"Identifier to uniquely specify a patient.",
        b"datatype": b"STRING",
        b"priority": b"1",
    }
    # Declared NUMBER: numeric. Declared STRING: text, though it holds "7".
    assert table.column("OS_MONTHS").to_pylist() == [12.5, 0.0]
    assert table.column("GLEASON").to_pylist() == ["7", None]
    assert table.column("GENDER").to_pylist() == ["Female", "Unknown"]


def test_declared_number_with_text_is_refused() -> None:
    text = CLINICAL_PATIENT.replace("P-2\tUnknown\t0", "P-2\tUnknown\tNA")
    with pytest.raises(ValueError, match="OS_MONTHS: declared NUMBER"):
        _read("data_clinical_patient.txt", text)


def test_cna_matrix_is_text() -> None:
    table = _read("data_cna.txt", CNA)
    assert {str(f.type) for f in table.schema} == {"string"}
    assert table.column("P-1-T02-IM7").to_pylist() == ["-2", "-1.5"]
    assert table.column("P-2-T01-IM6").to_pylist() == ["NA", "NA"]


def test_ragged_rows_are_refused() -> None:
    with pytest.raises(ValueError, match="line 3 has 1 cells"):
        _read("data_x.txt", "A\tB\n1\t2\n3\n")


def test_partial_clinical_header_is_refused() -> None:
    with pytest.raises(ValueError, match="2 '#' header lines"):
        _read("data_x.txt", "#a\tb\n#c\td\nA\tB\n1\t2\n")

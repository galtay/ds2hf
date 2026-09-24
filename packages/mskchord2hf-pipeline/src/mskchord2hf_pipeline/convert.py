"""One cBioPortal data file -> one Arrow table, values as published.

The license allows format changes and nothing more, so every rule here is
a statement about representation, never about meaning:

- **Lines.** A line ends at `\\n`, optionally preceded by `\\r`; the last
  line may lack one. Tabs separate cells; there is no quoting.
- **Missing.** An empty cell becomes null. Every other token (`NA`,
  `Unknown`, `None`, `-`) is a value and stays text, verbatim, padding
  spaces included.
- **Numbers.** A column becomes int64 when every non-empty cell is an
  integer that prints back to the same text, float64 when every non-empty
  cell is a decimal number float64 holds without losing a digit. Otherwise
  text. Clinical files declare a type per column; there, `STRING` is text
  whatever the values look like, and `NUMBER` must pass the numeric test
  or the build refuses.
- **Header lines.** Clinical files open with four `#` lines (display name,
  description, datatype, priority). They move into each field's parquet
  metadata, one key per line.
- **The CNA matrix** is text throughout: its `NA` marks a gene the sample's
  panel did not cover, and typing one sample column as numbers and the next
  as text would misrepresent a uniform matrix.

Row order, column order and column names are the file's.
"""

from __future__ import annotations

import re
from decimal import Decimal
from pathlib import Path

import pyarrow as pa

# The clinical header lines, in file order, as field-metadata keys.
HEADER_KEYS = ("display_name", "description", "datatype", "priority")

# Configs whose value columns stay text; see the module docstring.
TEXT_ONLY = frozenset({"cna"})

_INT = re.compile(r"-?(0|[1-9]\d*)\Z")
# No leading zeros and no bare leading dot: "007" and ".5" stay text, so a
# zero-padded identifier can never lose its padding.
_DECIMAL = re.compile(r"-?(0|[1-9]\d*)(\.\d*)?([eE][-+]?\d+)?\Z")


def config_name(member: str) -> str | None:
    """`data_timeline_surgery.txt` -> `timeline_surgery`; None for non-data files.

    Only the `data_` prefix and the extension go; a dot inside the stem
    (`data_cna_hg19.seg`) becomes an underscore because config names are
    path segments.
    """
    if "/" in member or not member.startswith("data_"):
        return None
    stem = member.removeprefix("data_")
    return stem.removesuffix(".txt").replace(".", "_")


def split_lines(text: str) -> list[list[str]]:
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return [(line[:-1] if line.endswith("\r") else line).split("\t") for line in lines]


def is_int(v: str) -> bool:
    return bool(_INT.match(v)) and str(int(v)) == v


def is_exact_float(v: str) -> bool:
    # repr is the shortest text that parses back to the same float, so this
    # holds exactly when float64 keeps every digit the text wrote.
    return bool(_DECIMAL.match(v)) and Decimal(repr(float(v))) == Decimal(v)


def _numeric_type(values: list[str]) -> pa.DataType | None:
    present = [v for v in values if v != ""]
    if not present:
        return None
    if all(is_int(v) for v in present):
        return pa.int64()
    if all(is_exact_float(v) for v in present):
        return pa.float64()
    return None


def _array(values: list[str], typ: pa.DataType) -> pa.Array:
    if typ == pa.int64():
        return pa.array([None if v == "" else int(v) for v in values], typ)
    if typ == pa.float64():
        return pa.array([None if v == "" else float(v) for v in values], typ)
    return pa.array([None if v == "" else v for v in values], pa.string())


def read_table(member: str, data: bytes) -> pa.Table:
    """Parse one data file from the study tarball into an Arrow table."""
    config = config_name(member)
    if config is None:
        raise ValueError(f"{member}: not a cBioPortal data file")
    rows = split_lines(data.decode("utf-8"))

    n_meta = 0
    while n_meta < len(rows) and rows[n_meta][0].startswith("#"):
        n_meta += 1
    if n_meta not in (0, len(HEADER_KEYS)):
        raise ValueError(f"{member}: {n_meta} '#' header lines; expected 0 or 4")
    meta = [[row[0][1:], *row[1:]] for row in rows[:n_meta]]
    header, body = rows[n_meta], rows[n_meta + 1 :]

    if len(set(header)) != len(header):
        raise ValueError(f"{member}: duplicate column names")
    for i, row in enumerate(body):
        if len(row) != len(header):
            n = n_meta + 2 + i
            raise ValueError(f"{member}: line {n} has {len(row)} cells, not {len(header)}")
    for line in meta:
        if len(line) != len(header):
            raise ValueError(f"{member}: header line has {len(line)} cells, not {len(header)}")

    fields, arrays = [], []
    for j, name in enumerate(header):
        values = [row[j] for row in body]
        field_meta: dict[str, str] | None = None
        if meta:
            field_meta = {key: line[j] for key, line in zip(HEADER_KEYS, meta, strict=True)}
            declared = field_meta["datatype"]
            if declared == "NUMBER":
                typ = _numeric_type(values)
                if typ is None:
                    raise ValueError(f"{member}.{name}: declared NUMBER but holds non-numbers")
            elif declared in ("STRING", "BOOLEAN"):
                typ = pa.string()
            else:
                raise ValueError(f"{member}.{name}: unknown declared datatype {declared!r}")
        elif config in TEXT_ONLY:
            typ = pa.string()
        else:
            typ = _numeric_type(values) or pa.string()
        fields.append(pa.field(name, typ, metadata=field_meta))
        arrays.append(_array(values, typ))
    schema = pa.schema(fields, metadata={"source_file": member})
    return pa.Table.from_arrays(arrays, schema=schema)


def write_table(table: pa.Table, path: Path) -> None:
    import pyarrow.parquet as pq

    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, path, compression="zstd")

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

REPO_ID = "gabrielaltay/msk-chord-2024"

# Sample-level configs and the column in each that names the sample. Every
# other sample-keyed config is nested under `samples` by these keys.
SAMPLE_CHILDREN: dict[str, str] = {
    "gene_panel_matrix": "SAMPLE_ID",
    "mutations": "Tumor_Sample_Barcode",
    "sv": "Sample_Id",
    "cna_hg19_seg": "ID",
}

# The CNA matrix's marker for a gene the sample's panel did not cover.
CNA_NOT_PROFILED = "NA"


def cna_long(wide: pa.Table, *, keep_unprofiled: bool = False) -> pa.Table:
    """The gene x sample CNA matrix as rows of (Hugo_Symbol, SAMPLE_ID, cna).

    `cna` is float64, since the matrix is not integer-only. Cells marked
    `NA`, genes the sample's panel did not cover, are dropped unless
    `keep_unprofiled`, which keeps them as null. Any other non-number
    raises rather than being guessed at.
    """
    genes = wide.column("Hugo_Symbol").combine_chunks()
    symbols, samples, values = [], [], []
    for name in wide.column_names[1:]:
        text = wide.column(name).combine_chunks()
        unprofiled = pc.fill_null(pc.equal(text, CNA_NOT_PROFILED), True)
        numbers = pc.cast(pc.if_else(unprofiled, None, text), pa.float64())
        if keep_unprofiled:
            symbols.append(genes)
            values.append(numbers)
        else:
            symbols.append(pc.filter(genes, pc.invert(unprofiled)))
            values.append(pc.filter(numbers, pc.invert(unprofiled)))
        samples.append(pa.repeat(name, len(values[-1])))
    return pa.table(
        {
            "Hugo_Symbol": pa.chunked_array(symbols, pa.string()),
            "SAMPLE_ID": pa.chunked_array(samples, pa.string()),
            "cna": pa.chunked_array(values, pa.float64()),
        }
    )


def _nest(
    parent_keys: pa.ChunkedArray, child: pa.Table, key: str, order_by: str | None = None
) -> pa.ListArray:
    """`child` rows grouped under `parent_keys`, one list per parent, in parent order.

    Within a list, rows keep their source order, or are sorted by `order_by`
    (nulls last) with source order breaking ties. A child key with no
    parent raises: silently dropping rows would misstate the data.
    """
    position = {k: i for i, k in enumerate(parent_keys.to_pylist())}
    keys = child.column(key).to_pylist()
    orphans = sorted({k for k in keys if k not in position}, key=str)
    if orphans:
        raise ValueError(f"{len(orphans)} {key} value(s) have no parent row, e.g. {orphans[:3]}")
    pos = [position[k] for k in keys]
    sort = pa.table(
        {"pos": pa.array(pos, pa.int64()), "row": pa.array(range(len(pos)), pa.int64())}
    )
    sort_keys = [("pos", "ascending"), ("row", "ascending")]
    if order_by is not None:
        sort = sort.append_column("by", child.column(order_by))
        sort_keys.insert(1, ("by", "ascending"))
    ordered = child.take(pc.sort_indices(sort, sort_keys=sort_keys, null_placement="at_end"))

    counts = [0] * len(position)
    for p in pos:
        counts[p] += 1
    offsets = [0]
    for c in counts:
        offsets.append(offsets[-1] + c)
    values = pa.StructArray.from_arrays(
        [col.combine_chunks() for col in ordered.columns], fields=list(ordered.schema)
    )
    return pa.ListArray.from_arrays(pa.array(offsets, pa.int32()), values)


class Chord:
    """The MSK-CHORD mirror: its configs as Arrow tables, plus local restructurings."""

    def __init__(self, root: Path | str):
        self.root = Path(root)
        self._tables: dict[str, pa.Table] = {}
        self._patients: dict[bool, pa.Table] = {}

    @classmethod
    def from_hub(cls, repo_id: str = REPO_ID, revision: str | None = None) -> Chord:
        """Download every config's parquet (not the source tarball) and open them."""
        from huggingface_hub import snapshot_download

        path = snapshot_download(
            repo_id, repo_type="dataset", revision=revision, allow_patterns=["*/data.parquet"]
        )
        return cls(path)

    @property
    def configs(self) -> list[str]:
        return sorted(p.parent.name for p in self.root.glob("*/data.parquet"))

    def table(self, config: str) -> pa.Table:
        """One config, exactly as published."""
        if config not in self._tables:
            self._tables[config] = pq.read_table(self.root / config / "data.parquet")
        return self._tables[config]

    def dictionary(self, config: str) -> list[dict[str, str]]:
        """Per-column display name, description, datatype and priority.

        Clinical configs only: these are the four `#` header lines of the
        cBioPortal clinical file, kept as parquet field metadata.
        """
        schema = pq.read_schema(self.root / config / "data.parquet")
        return [
            {"column": f.name, **{k.decode(): v.decode() for k, v in (f.metadata or {}).items()}}
            for f in schema
        ]

    @property
    def timelines(self) -> list[str]:
        return [c for c in self.configs if c.startswith("timeline_")]

    def patient_table(self, *, include_cna: bool = False) -> pa.Table:
        """One row per patient: clinical columns, `samples`, one list per timeline.

        - `samples` holds the patient's `clinical_sample` rows, each with a
          list per sample-level config (`gene_panel_matrix`, `mutations`,
          `sv`, `cna_hg19_seg`), and `cna` from `cna_long` if `include_cna`.
        - Each `timeline_*` column lists that timeline's events, sorted by
          `START_DATE`.

        Built locally; the license does not let you share the result.
        """
        if include_cna in self._patients:
            return self._patients[include_cna]
        samples = self.table("clinical_sample")
        sample_ids = samples.column("SAMPLE_ID")
        for config, key in SAMPLE_CHILDREN.items():
            samples = samples.append_column(config, _nest(sample_ids, self.table(config), key))
        if include_cna:
            samples = samples.append_column(
                "cna", _nest(sample_ids, cna_long(self.table("cna")), "SAMPLE_ID")
            )

        patients = self.table("clinical_patient")
        patient_ids = patients.column("PATIENT_ID")
        patients = patients.append_column("samples", _nest(patient_ids, samples, "PATIENT_ID"))
        for config in self.timelines:
            nested = _nest(patient_ids, self.table(config), "PATIENT_ID", order_by="START_DATE")
            patients = patients.append_column(config, nested)
        self._patients[include_cna] = patients
        return patients

    def patient(self, patient_id: str, *, include_cna: bool = False) -> dict[str, Any]:
        """One patient's row of `patient_table`, as a dict."""
        table = self.patient_table(include_cna=include_cna)
        rows = table.filter(pc.equal(table.column("PATIENT_ID"), patient_id)).to_pylist()
        if not rows:
            raise KeyError(patient_id)
        return rows[0]

    def patients(self, *, include_cna: bool = False, batch_size: int = 256) -> Iterator[dict]:
        """Every row of `patient_table`, as dicts, in `clinical_patient` order."""
        for batch in self.patient_table(include_cna=include_cna).to_batches(batch_size):
            yield from batch.to_pylist()

"""Generator for the clinical *_FIELDS lists in ds2hf/tcga/schema.py.

Reads each entity YAML from the gdcdictionary repo and prints a pyarrow field
list to stdout. Intended use: regenerate when the dictionary is updated, then
manually splice the output into `src/ds2hf/tcga/schema.py`.

Usage:
    uv run --with pyyaml python scripts/regenerate_clinical_fields.py \\
        /path/to/gdcdictionary > /tmp/clinical_fields.py

Then diff /tmp/clinical_fields.py against the corresponding sections of
src/ds2hf/tcga/schema.py and update.

Per-fetch dictionary snapshots at <data-dir>/raw/gdc_dictionary.<X>.<Y>.json
are the canonical record of what the GDC was serving when the parquet data
shipped — if they diverge from the static schema, regenerate.
"""

from __future__ import annotations

import sys
from pathlib import Path

import yaml

# System properties + relationship pointers we never want in our patient row.
SKIP_NAMES = {
    "type",
    "id",
    "state",
    "created_datetime",
    "updated_datetime",
    "project_id",
    "cases",
    "projects",
    "diagnoses",
    "samples",
    "portions",
    "analytes",
    "aliquots",
    "follow_ups",
    "treatments",
    "exposures",
    "family_histories",
    "demographic",
    "tissue_source_sites",
    "annotations",
    "molecular_tests",
    "other_clinical_attributes",
    "pathology_details",
    "slides",
    "read_groups",
}

# Entity dependency order: leaf first so parents can reference children.
DEPENDENCY_ORDER = [
    ("demographic", "DEMOGRAPHIC_FIELDS", []),
    ("treatment", "TREATMENT_FIELDS", []),
    ("diagnosis", "DIAGNOSIS_FIELDS", [("treatments", "TREATMENT_FIELDS")]),
    ("follow_up", "FOLLOW_UP_FIELDS", []),
    ("exposure", "EXPOSURE_FIELDS", []),
    ("family_history", "FAMILY_HISTORY_FIELDS", []),
    ("aliquot", "ALIQUOT_FIELDS", []),
    ("analyte", "ANALYTE_FIELDS", [("aliquots", "ALIQUOT_FIELDS")]),
    ("portion", "PORTION_FIELDS", [("analytes", "ANALYTE_FIELDS")]),
    ("sample", "SAMPLE_FIELDS", [("portions", "PORTION_FIELDS")]),
]

# Each entity's id field — system properties in GDC's dict but we always keep
# them as our schema's primary keys.
ID_FIELD = {ent: f"{ent}_id" for ent, _, _ in DEPENDENCY_ORDER}


def pa_type_for(name: str, prop_def: dict) -> str:
    """Map one dictionary property to a pyarrow type literal."""
    t = prop_def.get("type")
    if t is None and "oneOf" in prop_def:
        for opt in prop_def["oneOf"]:
            if isinstance(opt, dict) and "type" in opt:
                t = opt["type"]
                break
    if t is None and "enum" in prop_def:
        return "pa.string()"
    if isinstance(t, list):
        non_null = [x for x in t if x != "null"]
        t = non_null[0] if non_null else "string"
    if t == "integer":
        return "pa.int64()"
    if t == "number":
        return "pa.float64()"
    if t == "boolean":
        return "pa.bool_()"
    if t == "array":
        return "pa.list_(pa.string())"
    if t == "string":
        return "pa.string()"
    # Heuristics for refs that don't carry inline type info
    if "age_at" in name or name in {
        "year_of_birth",
        "year_of_diagnosis",
        "year_of_death",
        "year_of_follow_up",
    }:
        return "pa.int64()"
    return "pa.string()"


def main(dict_repo: Path) -> None:
    schemas_dir = dict_repo / "src/gdcdictionary/schemas"
    if not schemas_dir.is_dir():
        sys.exit(f"expected schemas dir at {schemas_dir}, not found")

    for ent, varname, children in DEPENDENCY_ORDER:
        s = yaml.safe_load((schemas_dir / f"{ent}.yaml").read_text())
        props = s.get("properties", {})
        sysprops = set(s.get("systemProperties", []))
        skip = SKIP_NAMES | (sysprops - {ID_FIELD[ent], "submitter_id"})

        fields = [(ID_FIELD[ent], "pa.string()"), ("submitter_id", "pa.string()")]
        for name in sorted(props):
            if name in skip or name in {ID_FIELD[ent], "submitter_id"}:
                continue
            defn = props[name]
            if not isinstance(defn, dict):
                continue
            fields.append((name, pa_type_for(name, defn)))

        print(f"{varname}: list[pa.Field] = [")
        for name, ptype in fields:
            print(f'    pa.field("{name}", {ptype}),')
        for child_name, child_var in children:
            print(f'    pa.field("{child_name}", pa.list_(pa.struct({child_var}))),')
        print("]\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: regenerate_clinical_fields.py <path-to-gdcdictionary-repo>")
    main(Path(sys.argv[1]))

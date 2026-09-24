"""Typed pydantic models + pyarrow schemas for the
`gabrielaltay/tcga-{patients,tabular}-open` HuggingFace datasets.

The pyarrow `*_FIELDS` lists in `ds2hf.tcga.schema` are the single source of
truth for the dataset shape (regenerated from gdcdictionary YAMLs); the
pydantic models in `ds2hf.tcga.models` are derived from those same lists, so
the two stay in sync by construction. See `TcgaHfPatient` for the typed
patient row + helper joins (tumor/normal pairs, mutations-by-gene,
expression lookup, longitudinal timeline).
"""

from ds2hf.tcga.models import TcgaHfPatient

__all__ = ["TcgaHfPatient"]

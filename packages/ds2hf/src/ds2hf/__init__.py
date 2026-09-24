"""Read-side companions for the ds2hf HuggingFace datasets, one module per source.

- `ds2hf.tcga`: typed pydantic models and pyarrow schemas for the TCGA datasets.
- `ds2hf.mskchord`: a local loader that builds per-patient views of the
  MSK-CHORD mirror on your machine.
"""

__version__ = "0.1.0"

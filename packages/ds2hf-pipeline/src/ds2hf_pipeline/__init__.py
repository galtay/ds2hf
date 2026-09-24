"""Build pipelines for the ds2hf HuggingFace datasets, one subpackage per source.

- `ds2hf_pipeline.tcga`: TCGA from the NCI Genomic Data Commons.
- `ds2hf_pipeline.mskchord`: MSK-CHORD from cBioPortal, format conversion only.

Sources share `hf_upload` and nothing else: each has its own licence, keys
and verification, and the code says so rather than abstracting over it.
"""

__version__ = "0.1.0"

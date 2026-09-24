"""Build pipeline that downloads public TCGA data from the NCI GDC and
publishes it as the `gabrielaltay/tcga-{patients,tabular}-open` HuggingFace
datasets, and their companions. The user-facing layer is
`ds2hf-pipeline tcga`; the read-side typed models and pyarrow schemas are
`ds2hf.tcga`.
"""

"""MSigDB gene-set collections (GMT), fetched and md5-pinned.

Companion to `ssgsea.py`, which consumes these sets. Same handling as the
Liu CDR workbook in `cdr.py`: download once into `<data-dir>/raw/msigdb/`,
verify against the hash we tested with, and surface a mismatch rather than
silently consuming new bytes.

Pinning matters more here than usual. MSigDB re-releases regularly (the
collections are curated, so set membership genuinely changes between
versions), and gene-set membership feeds straight into every score. A
dataset that says "ssGSEA of Hallmark" without naming the release is not
reproducible.

Only collections we can redistribute scores from are listed. MSigDB is
CC BY 4.0, but some constituent collections carry extra restrictions —
notably the KEGG-derived sets — which are therefore absent.

Which of the listed collections are actually scored is `SSGSEA_COLLECTIONS`
in `tcga2hf.schema`, not this dict: WikiPathways is fetchable and pinned
here but deliberately not scored (see its entry below).

Reference: https://www.gsea-msigdb.org/gsea/msigdb
License terms: https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import httpx

MSIGDB_VERSION = "2026.1.Hs"
_BASE_URL = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release"


@dataclass(frozen=True)
class Collection:
    """One MSigDB collection we can fetch, pinned by md5."""

    key: str
    stem: str
    md5: str
    description: str

    @property
    def file_name(self) -> str:
        return f"{self.stem}.v{MSIGDB_VERSION}.symbols.gmt"

    @property
    def url(self) -> str:
        return f"{_BASE_URL}/{MSIGDB_VERSION}/{self.file_name}"


# md5s are of the 2026-08 snapshots we validated against. The scored set
# was chosen by measuring gene coverage against our expression universe
# across 14 candidate collections rather than by reputation — all five
# exceed a 0.97 mean match fraction, while well-curated collections like
# C3:TFT (0.77) and C1 positional (0.39) were rejected because their genes
# largely fall outside a coding-gene matrix. See `dev_todo/ssGSEA.md`.
COLLECTIONS: dict[str, Collection] = {
    "hallmark": Collection(
        key="hallmark",
        stem="h.all",
        md5="367eec875967c2cfbf664a1a065b7b8d",
        description="MSigDB Hallmark — 50 coherent, deliberately non-redundant signatures",
    ),
    "reactome": Collection(
        key="reactome",
        stem="c2.cp.reactome",
        md5="1516b5d15611415d1996c92b7cb6d1cc",
        description="MSigDB C2:CP:REACTOME — 1,839 canonical pathways",
    ),
    "pid": Collection(
        key="pid",
        stem="c2.cp.pid",
        md5="291508046f73d82d13e5efb47492fa47",
        description="MSigDB C2:CP:PID — 196 NCI-Nature cancer signalling pathways",
    ),
    "oncogenic": Collection(
        key="oncogenic",
        stem="c6.all",
        md5="aba0e2214ff63327ae3fb0ce4bcd11c2",
        description=(
            "MSigDB C6 — 189 oncogenic signatures (oncogene / tumour-suppressor perturbation)"
        ),
    ),
    "cancer_cell_atlas": Collection(
        key="cancer_cell_atlas",
        stem="c4.3ca",
        md5="ff9902288655ff2ab88fcb5cbc4a95dd",
        description=(
            "MSigDB C4:3CA — 148 Curated Cancer Cell Atlas meta-programs (single-cell derived)"
        ),
    ),
    # Fetchable and pinned, but deliberately not in SSGSEA_COLLECTIONS: its
    # weakest sets are miRNA-centric and TCGA's poly-A RNA-Seq leaves miRNAs
    # 94.8% exactly zero, so 49 of its sets would score something other than
    # what their names claim. Revisit once the separate miRNA-Seq assay is in
    # the dataset, and re-measure rather than assume.
    "wikipathways": Collection(
        key="wikipathways",
        stem="c2.cp.wikipathways",
        md5="8e8a38972816a3997a557d6dd625138a",
        description="MSigDB C2:CP:WIKIPATHWAYS — 925 community-curated pathways",
    ),
}

# Every MSigDB gene set has a definition page at a URL derivable from its
# name — verified across all 26,937 pathways in the collections above. We
# template it onto each row rather than parsing the GMT's second field, the
# same way `gdc_portal_url` is templated from `case_id`, so the HF viewer
# renders a clickable link to the authoritative definition.
GENESET_URL_PREFIX = "https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/"


def geneset_url(pathway: str) -> str:
    return f"{GENESET_URL_PREFIX}{pathway}"


def md5_of(path: Path, chunk: int = 1 << 16) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def fetch_collection(
    collection: str,
    raw_dir: Path,
    *,
    verify: bool = True,
    timeout: float = 120.0,
) -> Path:
    """Download one collection's GMT into `<raw_dir>/msigdb/`; return the path.

    Re-runs are cheap: an existing file whose md5 matches the pin is left
    alone. With `verify=False` the hash is reported rather than enforced,
    which is how a new collection gets its pin recorded the first time.
    """
    if collection not in COLLECTIONS:
        raise ValueError(f"unknown collection {collection!r}; known: {sorted(COLLECTIONS)}")
    spec = COLLECTIONS[collection]
    out_dir = raw_dir / "msigdb"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / spec.file_name

    if out_path.exists():
        digest = md5_of(out_path)
        if not verify or not spec.md5 or digest == spec.md5:
            return out_path
        raise ValueError(
            f"{out_path} md5 {digest} != pinned {spec.md5}. MSigDB may have re-released "
            f"{MSIGDB_VERSION}; inspect before overwriting the pin."
        )

    with httpx.Client(timeout=timeout, follow_redirects=True) as client:
        resp = client.get(spec.url)
        resp.raise_for_status()
        out_path.write_bytes(resp.content)

    digest = md5_of(out_path)
    if verify and spec.md5 and digest != spec.md5:
        out_path.unlink()
        raise ValueError(
            f"downloaded {spec.file_name} has md5 {digest}, expected {spec.md5}. "
            "Refusing to keep unverified gene sets."
        )
    return out_path

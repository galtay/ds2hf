from __future__ import annotations

import re
from pathlib import Path

import pytest

CARD_MODULE = Path(__file__).resolve().parents[1] / "src" / "tcga2hf_pipeline" / "dataset_card.py"

# Markdown link targets and bare URLs in the card templates. Trailing
# punctuation that belongs to the prose rather than the URL is stripped.
_URL = re.compile(r"https?://[^\s\)\]\"'<>]+")


def card_urls() -> list[str]:
    """Every fixed external URL the generated dataset cards can contain.

    URLs holding an f-string placeholder — `.../datasets/{repo}/resolve/...`
    — are skipped: they are only a URL once a card is rendered, and what
    they resolve to depends on the repo being written.
    """
    found = _URL.findall(CARD_MODULE.read_text())
    return sorted({u.rstrip(".,;:") for u in found if "{" not in u})


def test_cards_cite_the_gdc_expression_pipeline() -> None:
    """The strandedness guidance rests on GDC's own statement, so the quote
    and the link backing it both have to survive edits to the card.

    Matched against whitespace-collapsed source rather than the literal
    text: an earlier version asserted on the exact line break and broke the
    moment the prose was rewrapped, which says nothing about the citation.
    """
    text = CARD_MODULE.read_text()
    collapsed = " ".join(text.replace(">", " ").split())
    assert "[gdc-mrna]: https://docs.gdc.cancer.gov" in text
    assert "all RNA-Seq reads are treated as unstranded during analyses" in collapsed
    assert "mRNA Analysis Pipeline][gdc-mrna], Introduction" in collapsed


def test_every_referenced_link_is_defined() -> None:
    """A `[text][label]` with no `[label]: url` renders as literal brackets.

    Cheap to get wrong when editing prose, and invisible until someone reads
    the published card.
    """
    text = CARD_MODULE.read_text()
    defined = set(re.findall(r"^\[([a-z0-9-]+)\]: \S+", text, re.MULTILINE))
    used = set(re.findall(r"\]\[([a-z0-9-]+)\]", text))
    assert used - defined == set(), f"undefined link labels: {sorted(used - defined)}"


@pytest.mark.network
def test_card_urls_resolve() -> None:
    """Every URL a card publishes must actually load.

    A 404 in a dataset card is a slow failure: it ships to the Hub and stays
    wrong until a reader clicks it. `docs.gdc.cancer.gov` has reorganised
    Encyclopedia paths before -- `/pages/Biospecimen/` became
    `/pages/Biospecimen_Data/` -- and the broken link sat on every published
    card until this test existed.
    """
    import httpx

    broken: list[str] = []
    for url in card_urls():
        try:
            response = httpx.get(url, timeout=30.0, follow_redirects=True)
        except httpx.HTTPError as exc:  # network flake, not a broken card
            pytest.skip(f"could not reach {url}: {type(exc).__name__}")
        if response.status_code >= 400:
            broken.append(f"{response.status_code} {url}")
    assert not broken, "dead links in dataset cards:\n  " + "\n  ".join(broken)


def test_gdc_references_defaults_to_every_doc() -> None:
    from tcga2hf_pipeline.dataset_card import _GDC_DOCS, _gdc_references

    section = _gdc_references()
    assert section.count("\n- ") == len(_GDC_DOCS)


def test_gdc_references_keeps_the_caller_s_order() -> None:
    from tcga2hf_pipeline.dataset_card import _gdc_references

    section = _gdc_references("barcode", "maf")
    assert section.index("Barcode") < section.index("MAF")
    assert "Biospecimen" not in section


def test_gdc_references_rejects_an_unknown_key() -> None:
    """A typo'd key would otherwise drop a reference silently."""
    from tcga2hf_pipeline.dataset_card import _gdc_references

    with pytest.raises(KeyError, match="unknown GDC doc key"):
        _gdc_references("mrna", "maff")


def test_expression_card_omits_irrelevant_specs() -> None:
    """The expression dataset carries no MAF and no biospecimen tree."""
    import re as _re

    text = CARD_MODULE.read_text()
    call = _re.search(r"\+ _gdc_references\(([^)]*)\)", text)
    assert call, "expression card no longer calls _gdc_references with keys"
    keys = call.group(1)
    assert "mrna" in keys
    assert "maf" not in keys
    assert "biospecimen" not in keys


def test_frontmatter_is_one_key_per_line() -> None:
    """YAML frontmatter must not be collapsed into a single line.

    HF reads `configs:` out of this block to decide what to serve, so a
    malformed one takes the dataset viewer down entirely. A prose-unwrapping
    pass once joined every key onto one line here, and a whitespace-
    normalising comparison reported the cards unchanged, because collapsing
    whitespace is exactly what hides this.
    """
    text = CARD_MODULE.read_text()
    blocks = re.findall(r"^license: other.*$", text, re.MULTILINE)
    assert blocks, "no frontmatter templates found"
    for block in blocks:
        assert block == "license: other", f"frontmatter collapsed: {block[:80]}"


@pytest.mark.parametrize("key", ["license_name", "license_link", "pretty_name", "tags"])
def test_frontmatter_keys_start_their_own_line(key: str) -> None:
    text = CARD_MODULE.read_text()
    assert re.search(rf"^{key}:", text, re.MULTILINE), f"{key} does not begin a line"

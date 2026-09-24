"""Local loader for the MSK-CHORD HuggingFace mirror.

The mirror publishes each cBioPortal file as its own config, values as
published: the study's CC BY-NC-ND 4.0 license lets it be shared only
unmodified. This package restructures it on your machine instead, into a
long-format copy-number table and one nested row per patient.

**Do not share what these functions return.** The license (section
2(a)(1)(B)) lets you produce restructured versions for non-commercial use,
but not share them. Share this code, or the unmodified mirror, instead.

Usage:

    from mskchord2hf import Chord

    chord = Chord.from_hub()
    patients = chord.patient_table()          # one row per patient, nested
    p = chord.patient("P-0000012")            # the same row as a dict
    for sample in p["samples"]:
        print(sample["SAMPLE_ID"], [m["Hugo_Symbol"] for m in sample["mutations"]])
"""

from mskchord2hf.loader import REPO_ID, Chord, cna_long

__all__ = ["REPO_ID", "Chord", "cna_long"]

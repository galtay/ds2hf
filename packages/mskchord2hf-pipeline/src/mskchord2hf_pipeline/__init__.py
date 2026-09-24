"""Mirror the MSK-CHORD cBioPortal study onto the HuggingFace Hub.

MSK-CHORD is licensed CC BY-NC-ND 4.0: it may be shared unmodified, and the
license's section 2(a)(4) exempts format changes from counting as
modification. So this package converts formats and does nothing else: one
config per source file, every value and row as published. Anything that
restructures the data (joins, per-patient nesting, reshaping) lives in the
`mskchord2hf` loader and runs on the user's own machine.
"""

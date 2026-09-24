"""A miniature msk_chord_2024 study: every file shape the real one has."""

from __future__ import annotations

import io
import tarfile
from pathlib import Path

from ds2hf_pipeline.mskchord import source

CLINICAL_PATIENT = (
    "#Patient Identifier\tSex\tOverall Survival (Months)\tGleason\n"
    "#Identifier to uniquely specify a patient.\tGender at birth\tOS from sequencing\tNLP score\n"
    "#STRING\tSTRING\tNUMBER\tSTRING\n"
    "#1\t970\t0\t750\n"
    "PATIENT_ID\tGENDER\tOS_MONTHS\tGLEASON\n"
    "P-1\tFemale\t12.5\t7\n"
    "P-2\tUnknown\t0\t\n"
)
CLINICAL_SAMPLE = (
    "#Sample Identifier\tPatient Identifier\tCancer Type\n"
    "#Sample id\tPatient id\tOncotree type\n"
    "#STRING\tSTRING\tSTRING\n"
    "#1\t1\t3000\n"
    "SAMPLE_ID\tPATIENT_ID\tCANCER_TYPE\n"
    "P-1-T01-IM6\tP-1\tBreast Cancer\n"
    "P-1-T02-IM7\tP-1\tNA\n"
    "P-2-T01-IM6\tP-2\tLung Cancer\n"
)
# CRLF line endings, padded units and no final newline, as the lab timelines ship.
TIMELINE_LABS = (
    "PATIENT_ID\tSTART_DATE\tSTOP_DATE\tEVENT_TYPE\tRESULT\tLR_UNIT_MEASURE\r\n"
    "P-1\t30\t\tLab_test\t2\tng/ml     \r\n"
    "P-1\t-5\t\tLab_test\t1.25\tng/ml     \r\n"
    "P-2\t0\t\tLab_test\t0\tng/ml     "
)
MUTATIONS = (
    "Hugo_Symbol\tTumor_Sample_Barcode\tStart_Position\tChromosome\tReference_Allele\n"
    "TP53\tP-1-T01-IM6\t7577120\t17\t-\n"
    "KRAS\tP-2-T01-IM6\t25398284\tX\tC\n"
)
SV = "Sample_Id\tSite1_Hugo_Symbol\tSV_Length\nP-1-T02-IM7\tALK\t\n"
SEG = "ID\tchrom\tloc.start\tseg.mean\nP-1-T01-IM6\t1\t100\t-4e-04\nP-2-T01-IM6\tX\t200\t0.25\n"
PANELS = (
    "SAMPLE_ID\tmutations\tcna\n"
    "P-1-T01-IM6\tIMPACT468\tIMPACT468\n"
    "P-1-T02-IM7\tIMPACT505\tIMPACT505\n"
    "P-2-T01-IM6\tIMPACT468\tNA\n"
)
CNA = "Hugo_Symbol\tP-1-T01-IM6\tP-1-T02-IM7\tP-2-T01-IM6\nTP53\t0\t-2\tNA\nMYC\t2\t-1.5\tNA\n"

MEMBERS = {
    "LICENSE": "Attribution-NonCommercial-NoDerivatives 4.0 International\n",
    "meta_study.txt": (
        "cancer_study_identifier: msk_chord_2024\nname: MSK-CHORD (MSK, Nature 2024)\n"
    ),
    "meta_cna.txt": "profile_description: Values: -2 = deep loss; 0 = neutral; 2 = amp.\n",
    "case_lists/cases_all.txt": "stable_id: msk_chord_2024_all\n",
    "data_clinical_patient.txt": CLINICAL_PATIENT,
    "data_clinical_sample.txt": CLINICAL_SAMPLE,
    "data_timeline_labs.txt": TIMELINE_LABS,
    "data_mutations.txt": MUTATIONS,
    "data_sv.txt": SV,
    "data_cna_hg19.seg": SEG,
    "data_gene_panel_matrix.txt": PANELS,
    "data_cna.txt": CNA,
}


def write_tarball(path: Path, members: dict[str, str]) -> Path:
    with tarfile.open(path, "w:gz") as tar:
        for name, text in members.items():
            data = text.encode()
            info = tarfile.TarInfo(f"{source.STUDY_ID}/{name}")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return path

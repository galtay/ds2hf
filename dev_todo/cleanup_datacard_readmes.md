Lets cleanup the HF dataset cards.

Goals:

* keep it minimal and scoped to things that will not change (e.g. link to the top level repo instead of specific files as the links are already stale).
* describe what we did and how it was constructed. leave usage instructions minimal


Things to keep,

* a core component that is shared between the patient and tabular datasets
* GDC References

Thinks to update,

* we dont need explicit lists of all TCGA projects included (twice). the assumption is that we have done this for all TCGA projects

* things like this are misleading "Source: GDC /cases endpoint, open-access tier only." we use much more than the cases endpoint and its a detail not worth a top level spot in the data card

* things like this that are more like things we discovered along the way don't need to be included. "Every days_to_* field anchors to the case's index_date (TCGA: almost always "Diagnosis"), per the dictionary, so clinical and biospecimen events share a single timeline.".

* this is a lot of text and combines two things "Lifted expression QC counts. Each Gene Expression Quantification record has the STAR per-feature quality-control counts N_unmapped, N_multimapping, N_noFeature, N_ambiguous lifted from the source Tab-Separated Values (TSV) file onto the row as scalar fields. The stranded_first / stranded_second columns are dropped — the GDC pipeline harmonizes by treating all RNA-Seq reads as unstranded, so unstranded is the canonical column." The first is trival and we don't have to mention it. The second does deserve its own bullet point

* the enrichment with rederived survival endpoints using the Liu 2018 methods is a big addition. It should have its own section but it should be after we describe the core TCGA data. we also don't need to go into the details of teh analysis, we can succintly describe what we did and point to the github repo for details.

* for "How this dataset is built" we don't need to go into the details of post calls or relist every TCGA project. However, we should describe the source of everything in each dataset and keep things like what filters were used with the data_type, data_format etc. Lets add the clinical supplements to that

* in the disclaimer we should mention that this project is not affiliated with GDC or NCI and is an experimental open source project that will potentially be unstable and change often

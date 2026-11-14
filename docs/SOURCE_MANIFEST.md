# Source selection and caveats

- `01_deseq2_validation.R`: historical standard DESeq2 verification, not threshold-aware archived p-value regeneration. Runs in its original validation bundle directory.
- `02_ortholog_venn_validation.R`: first-match Orthofinder gene-set validation; enrichment portion must run in the specified environment.
- `04_go_kegg_recovered.R`: historical compareCluster workflow; redundant second setReadable removed. **Do not use this file as the authoritative GO reproduction until it is tested under Bioconductor 3.21.**
- `05_collectri_validated_no_batch.R`: extracted from the tested Ensembl/no-batch branch, using `times=1000` and SYMBOL-to-ENSEMBL network targets. Needs clean-session validation.
- `scripts/historical_source/`: six supplied scripts preserved verbatim for provenance, including failed symbol-based and experimental batch variants. Do not execute as a pipeline.
- Historical S7/S9 GO term identities, ratios, geneID sets, counts and p-values were validated exactly against Bioconductor 3.21 output in the user's Mac run.
- Figures and supplemental methods are not yet consolidated.

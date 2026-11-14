# Manuscript Figure Reproducibility

Validated analytical components:

- DESeq2: 26,543 gene IDs per contrast, four contrasts; small numerical differences from archived standard DESeq2 results. Archived threshold-aware lists reproduce historical significance sets.
- Orthofinder/Venn: all eight human Ensembl gene sets match the historical workbooks exactly.
- GO: Bioconductor 3.21 default compareCluster reproduces all 630 historical condition-term records, including gene sets, ratios, counts, and p-values to machine precision.
- CollecTRI: Ensembl/no-batch workflow reproduces 724 TFs and S10/S12 statistics to machine precision.
- KEGG: 79 of 80 historical significant condition-pathway records retained; see KEGG_VERSION_NOTE.md.

Historical figure scripts were reviewed; some are exploratory, incomplete, or from other developmental stages. This repository reproduces audited analytical results but does not yet regenerate every final assembled manuscript figure. Publication figure files and source provenance should be archived separately after permission review.

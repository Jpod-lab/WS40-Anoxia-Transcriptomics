# KEGG Annotation Reproducibility

The archived threshold-aware gene sets and first-match human Orthofinder mappings were used for KEGG enrichment. Of 80 historically significant condition-pathway records, 79 remained significant with the tested current KEGG annotations.

The exception is antifolate resistance (hsa01523), 24hR upregulated:

| Measure | Historical | Current |
|---|---|---|
| Contributing genes | IL1B, ABCC2 | IL1B, ABCC2 |
| GeneRatio | 2/16 | 2/17 |
| BgRatio | 30/9447 | 30/9472 |
| BH-adjusted p-value | 0.0461 | 0.05314109 |

Both contributing genes remain associated with hsa01523. Changes in effective input coverage and annotation background are consistent with KEGG annotation drift. The exact historical KEGG snapshot and specific causal annotation change have not been identified. The archived result is retained as the manuscript reference; the current result is reported transparently.

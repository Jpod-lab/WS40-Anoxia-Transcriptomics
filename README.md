# WS40 12-dpd manuscript reproducibility
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24267894.svg)](https://doi.org/10.5281/zenodo.24267894)

**Status: repository staging package, not yet a complete executable reproduction.** This repository is intended to preserve the manuscript-era results and independently rerun them. Do not replace archived supplemental tables with modern annotation outputs.

## Validated checkpoints

- Original 12-dpd expression: 30 samples, five conditions with six replicates each.
- Historical threshold-aware DESeq2 up/down gene sets and Venn sets: exact agreement.
- Eight Venn-derived unique human Ensembl gene sets: exact agreement with first Orthofinder match per killifish locus.
- CollecTRI/decoupleR: 724 TFs, S10 means and S12 ANOVA statistics reproduced to machine precision using Ensembl expression IDs, network SYMBOL-to-ENSEMBL targets, centered VST, no inferred batch correction.
- GO: Bioconductor 3.21 default compareCluster(ont='ALL', keyType='ENSEMBL') reproduces all 451 Up and 179 Down historical condition-term IDs, with no extras. **Numerical column agreement still requires validation.**
- KEGG: 49/50 Up and 30/30 Down historical significant records reproduced; missing 24hR antifolate resistance hsa01523 remains present at current q=0.05314109 versus historical q=0.0461.

## Analysis stages

1. Original `se12dpd` and archived threshold-aware DESeq2 tables.
2. First-match Orthofinder mapping, human Ensembl gene sets, Venn validation.
3. GO via clusterProfiler compareCluster in Bioconductor 3.21 using historical workbook inputs, unsimplified table exports.
4. KEGG via `bitr(ENSEMBL -> ENTREZID)` then compareCluster.
5. GO Wang simplify(cutoff=0.5) only for designated figure outputs, not S7/S9.
6. CollecTRI/decoupleR with SYMBOL network targets converted to ENSEMBL; preserve validated no-batch results.

## Repository data policy

- Do not commit patient data, private work directories, proprietary resources, or raw sequencing files without confirming distribution rights.
- Add only approved original manuscript source tables, public reference data, and scripts.
- Store large raw inputs externally with checksums and retrieval instructions.
- Verify all manuscript numerical columns and figure panels before claiming full reproduction.

## Next required steps

- Copy actual validated source R scripts and historical workbooks into the appropriate folders (do not reconstruct code from memory).
- Run `scripts/validate_go_tables.R` on the original historical S7/S9 and Bioc 3.21 GO outputs.
- Confirm original manuscript figures against source outputs, especially Figure 5 and supplemental figures.
- Pin Docker image by digest and archive sessionInfo/package versions.
- Publish only after reviewing licensing and data sensitivity.

## Publication checklist
See `docs/RELEASE_CHECKLIST.md`. The Dockerfile pins the Bioconductor release tag but must be built and digest-pinned before publication. The numerical GO validation is confirmed by the user's R run (630/630 records, ratios and counts exact, p-values at machine precision); the contributing geneID field remains to be audited.

## Source-code consolidation
Six user-supplied scripts are preserved under `scripts/historical_source/`. Candidate cleaned scripts are staged under `scripts/` and require clean-run validation before being called production-ready. See `docs/SOURCE_MANIFEST.md`.

## Final run and figure guidance
See `docs/RUN_GUIDE.md`, `docs/FIGURE_VALIDATION_STATUS.md`, `docs/KEGG_VERSION_NOTE.md`, and `docs/FIGURE_REPRODUCIBILITY.md`. Inputs and historical tables are excluded from Git pending PI approval; the repository is not self-contained without them.

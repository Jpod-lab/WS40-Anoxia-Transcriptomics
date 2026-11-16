# WS40 12-dpd manuscript reproducibility

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24267894.svg)](https://doi.org/10.5281/zenodo.24267894)

**Status: Public manuscript reproducibility release (v1.0.1).** This repository provides portable analysis and manuscript figure-generation scripts, environment specifications, and validation documentation for the WS40 12-days-post-diapause transcriptomics study. Historical analytical results have been independently reproduced as described below. A fresh clone, Docker build, core package loading, and syntax checks passed. **Full end-to-end execution from a clean clone has not yet been demonstrated**, because historical input datasets and reference workbooks are distributed separately.

## Validated scientific checkpoints

- **Experimental design:** Original 12-dpd expression dataset contains 30 samples: five conditions with six biological replicates each.
- **DESeq2:** Reconstructed standard differential-expression analyses matched all 106,172 gene classifications across four contrasts (26,543 genes per contrast). Small numerical differences remain in some statistics; archived threshold-aware results are the reference for historical significant gene sets.
- **Ortholog mapping and Venn:** All eight historical up/down gene sets and the resulting unique human Ensembl sets were reproduced exactly using the first Orthofinder match per killifish locus.
- **GO enrichment:** With Bioconductor 3.21 and the validated `compareCluster` workflow, all 630 historical condition-term records (451 Up, 179 Down) were reproduced, including GeneRatio, BgRatio, Count, geneID, and p-values to machine precision. Historical S7/S9 outputs should not be replaced by results from newer annotation databases.
- **KEGG enrichment:** 79 of 80 historically significant records were recovered (49/50 Up; 30/30 Down). The Up 24hR antifolate-resistance pathway (`hsa01523`) has current q = 0.05314109 versus historical q = 0.0461. This is consistent with changed annotation/input mapping, but the precise cause has not been conclusively established. See `docs/KEGG_VERSION_NOTE.md`.
- **CollecTRI/decoupleR:** All 724 transcription factors and the historical S10/S12 results were reproduced to machine precision using Ensembl expression identifiers, SYMBOL-to-ENSEMBL network target mapping, and no inferred batch correction.
- **Figures:** Portable scripts for Venn, PCA, combined GO/KEGG, volcano, and transcription-factor figures are included. Selected regenerated figures were visually reviewed; output files are not claimed to be byte-identical to historical originals.

## Repository contents

| Location | Purpose |
| --- | --- |
| `scripts/01_deseq2_portable.R` | DESeq2 reconstruction and comparison |
| `scripts/02_ortholog_venn_portable.R` | Ortholog mapping and Venn set validation |
| `scripts/04_go_bioc321_portable.R` | Historical GO reproduction with Bioconductor 3.21 |
| `scripts/05_collectri_portable.R` | CollecTRI transcription-factor analysis |
| `scripts/figures/` | Manuscript figure-generation scripts |
| `scripts/validate_go_complete.R`, `scripts/validate_go_tables.R` | GO validation utilities |
| `environment/Dockerfile` | Bioconductor 3.21 environment definition |
| `docs/` | Run guide, provenance, figure status, and validation notes |
| `validation/` | Tracked validation summaries |

See `scripts/00_run_order.txt`, `docs/RUN_GUIDE.md`, `docs/SOURCE_MANIFEST.md`, `docs/FIGURE_REPRODUCIBILITY.md`, and `docs/FIGURE_VALIDATION_STATUS.md` for detailed guidance.

## Reproduction and environment

1. Clone the repository and review `docs/RUN_GUIDE.md`.
2. Obtain the historical inputs and supplemental workbooks through the approved data-availability process; stage them in the repository-relative paths expected by the scripts. They are intentionally not committed to Git.
3. Build the environment from the repository root:

   `docker build --platform linux/amd64 -t ws40-repro:bioc321 -f environment/Dockerfile .`

4. Follow `scripts/00_run_order.txt` and `docs/RUN_GUIDE.md` to execute the relevant analyses after inputs are staged.
5. Compare outputs with historical tables using the supplied validation utilities and documented acceptance criteria.

**Clean-environment checks completed:** independent Git clone, Docker image build, loading of DESeq2, clusterProfiler, org.Hs.eg.db and decoupleR, and successful R syntax parsing of all 11 tracked R scripts. These checks do not substitute for full execution with historical inputs. The Docker base image uses a release tag and has not been pinned by immutable digest; figure-specific dependencies and complete clean-run behavior have not all been validated.

## Data availability and safeguards

Historical input datasets, large raw sequencing files, and supplemental workbooks are excluded from this Git repository pending distribution approval. Do not commit restricted, proprietary, or otherwise non-distributable data. Preserve checksums and retrieval instructions for separately distributed inputs. Reproducing an analysis with updated annotation resources may change results, especially enrichment statistics; use the documented historical environment for manuscript-era comparisons.

## Citation and licensing

**Software author:** Patrick R. Clouser. Scientific supervision and project leadership were provided by Jason Podrabsky and the Podrabsky Lab. See `CITATION.cff` for machine-readable software citation metadata.

- **Archived software release:** v1.0.1, [Zenodo DOI: 10.5281/zenodo.24267894](https://doi.org/10.5281/zenodo.24267894)
- **Source repository:** [Jpod-lab/WS40-Anoxia-Transcriptomics](https://github.com/Jpod-lab/WS40-Anoxia-Transcriptomics)
- **License:** MIT; see `LICENSE`. The license covers the repository software, not automatically any separately obtained input datasets.

The archived v1.0.1 snapshot is fixed. Later documentation commits on `main` are not part of that archived snapshot; cite the archived version for exact software provenance.

## Remaining limitations

- Full end-to-end clean-clone analysis execution requires separately staged historical inputs and remains unverified.
- The Docker image has not been digest-pinned.
- Not all figure scripts and figure-specific dependencies have undergone clean-environment execution testing.
- The KEGG `hsa01523` historical/current significance difference remains documented rather than definitively resolved.
- Regenerated figures are not asserted to be byte-identical to historical files.

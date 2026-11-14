# Local reproduction guide

Run from repository root. Inputs under data/ and historical xlsx tables are deliberately gitignored and must be obtained through the approved data-availability process.

1. `WS40_REPO="$PWD" Rscript scripts/01_deseq2_portable.R`
2. `WS40_REPO="$PWD" Rscript scripts/02_ortholog_venn_portable.R`
3. GO uses Bioconductor 3.21 in Docker: `docker run --rm --platform linux/amd64 -v "$PWD:/work/repo" bioconductor/bioconductor_docker:RELEASE_3_21 Rscript /work/repo/scripts/04_go_bioc321_portable.R` (packages may be installed in disposable container).
4. `WS40_REPO="$PWD" Rscript scripts/05_collectri_portable.R`
5. Figure scripts under scripts/figures/ can be run with `WS40_REPO="$PWD" Rscript <script>` after installing their dependencies.

Use `scripts/validate_go_complete.R` with S7, Up GO, S9, Down GO arguments to audit GO. Historical GO 630/630 exact; CollecTRI S10/S12 machine precision. Standard DESeq2 has small numerical differences; threshold-aware archived results are authoritative for Venn sets. KEGG antifolate resistance hsa01523 remains a documented threshold-crossing difference.

The Dockerfile is not yet digest-pinned or clean-build validated. Do not claim that all figure files are byte-identical.

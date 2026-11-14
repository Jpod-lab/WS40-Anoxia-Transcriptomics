repo <- normalizePath(Sys.getenv("WS40_REPO", unset="."))
input_dir <- file.path(repo,"data/input/deseq2")
out_dir <- file.path(repo,"results/reproduced/deseq2")
# Run from this directory on a Mac with DESeq2 installed.
# This script intentionally uses the historical five-count total filter.
suppressPackageStartupMessages(library(DESeq2))
counts <- read.csv(file.path(input_dir,"WS40_12dpd_raw_counts.csv"), row.names=1, check.names=FALSE)
metadata <- read.csv(file.path(input_dir,"WS40_12dpd_sample_metadata.csv"), row.names=1, check.names=FALSE)
stopifnot(ncol(counts)==nrow(metadata), all(counts >= 0), all(counts == round(as.matrix(counts))))
# CSV export lost the original assay column names (V1..V30); restore sample order explicitly.
colnames(counts) <- metadata$Sample
rownames(metadata) <- metadata$Sample
stopifnot(identical(colnames(counts),rownames(metadata)))
metadata$treatment <- factor(metadata$treatment, levels=c("t0","SA","LA","SR","LR"))
stopifnot(!anyNA(metadata$treatment))
counts <- counts[rowSums(counts)>=5,,drop=FALSE]
dds <- DESeqDataSetFromMatrix(countData=round(as.matrix(counts)), colData=metadata, design=~treatment)
dds <- DESeq(dds)
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
for (condition in c("SA","LA","SR","LR")) {
  result <- results(dds, contrast=c("treatment",condition,"t0"), alpha=0.05, pAdjustMethod="BH")
  reference <- read.csv(file.path(input_dir,paste0("12dpd_",condition,"_t0_fullset.csv")),row.names=1,check.names=FALSE)
  result <- as.data.frame(result)
  common <- intersect(rownames(reference),rownames(result))
  fields <- intersect(c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"),names(reference))
  comparison <- data.frame(GeneID=common)
  for (field in fields) {
    comparison[[paste0(field,"_original")]] <- reference[common,field]
    comparison[[paste0(field,"_reproduced")]] <- result[common,field]
    comparison[[paste0(field,"_absdiff")]] <- abs(reference[common,field]-result[common,field])
  }
  write.csv(comparison,file.path(out_dir,paste0(condition,"_gene_by_gene_comparison.csv")),row.names=FALSE)
  write.csv(result,file.path(out_dir,paste0(condition,"_reproduced_fullset.csv")))
  cat(condition,"original",nrow(reference),"reproduced",nrow(result),"matched IDs",length(common),"\n")
  for (field in fields) cat(field,"max abs diff",max(comparison[[paste0(field,"_absdiff")]],na.rm=TRUE),"\n")
}
writeLines(capture.output(sessionInfo()),file.path(out_dir,"sessionInfo.txt"))

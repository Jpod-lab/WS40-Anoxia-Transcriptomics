FROM bioconductor/bioconductor_docker:RELEASE_3_21
RUN Rscript -e 'BiocManager::install(c("clusterProfiler","org.Hs.eg.db","GO.db","GOSemSim","AnnotationDbi","readxl","decoupleR","limma","DESeq2"),ask=FALSE,update=FALSE)'

# Final manuscript-aligned WS40 GO/KEGG rerun
# Exact archived threshold-aware lists -> first Orthofinder row per killifish locus
# -> unique human Ensembl IDs for GO -> bitr to unique Entrez IDs for KEGG
# -> unsimplified compareCluster tables; simplify only for GO figure candidates.
suppressPackageStartupMessages({library(clusterProfiler);library(org.Hs.eg.db);library(GOSemSim);library(readxl);library(dplyr)})
repo <- normalizePath(Sys.getenv('WS40_REPO', unset='.')); base <- file.path(repo,'data/input/historical_threshold')
out <- file.path(repo,'results/reproduced/ortholog_venn');dir.create(out,recursive=TRUE,showWarnings=FALSE)
labels <- c(SA='4hA',LA='24hA',SR='2hR',LR='24hR')
# Archived threshold-aware results are source of truth for manuscript gene lists.
paths <- setNames(vapply(names(labels),function(code)file.path(base,paste0('12dpd_',code,'_t0_res_sig_fullset.csv')),character(1)),names(labels))
stopifnot(all(file.exists(paths)))
res <- lapply(paths,function(f)read.csv(f,row.names=1,check.names=FALSE))
select_loci <- function(x,dir){good<-!is.na(x$padj)&x$padj<0.05&!is.na(x$log2FoldChange);good<-good & if(dir=='Up')x$log2FoldChange>0.5 else x$log2FoldChange< -0.5;unique(rownames(x)[good])}
loci <- lapply(c('Up','Down'),function(dir){z<-lapply(res,select_loci,dir=dir);names(z)<-labels[names(z)];z});names(loci)<-c('Up','Down')
# First Orthofinder record per killifish locus, preserving file order.
ortho_candidates<-file.path(repo,'data/input/Orthofinder_blast_merge_Hs_230628_MapIDS.csv')
op<-ortho_candidates[file.exists(ortho_candidates)][1];if(is.na(op))stop('Orthofinder file missing')
ortho<-read.csv(op,fileEncoding='latin1',stringsAsFactors=FALSE,check.names=FALSE)
first_map<-ortho[!duplicated(ortho$Locus_Al),]
first_map$EnsemblID_Hs<-sub('\\..*$','',as.character(first_map$EnsemblID_Hs))
map_ens<-function(ids){z<-first_map$EnsemblID_Hs[first_map$Locus_Al%in%ids];sort(unique(z[!is.na(z)&nzchar(z)&z!='NULL']))}
ens<-lapply(loci,function(s)lapply(s,map_ens))
map_entrez<-function(ids){if(!length(ids))return(character());z<-suppressMessages(bitr(ids,fromType='ENSEMBL',toType='ENTREZID',OrgDb=org.Hs.eg.db,drop=TRUE));sort(unique(as.character(z$ENTREZID[!is.na(z$ENTREZID)&nzchar(z$ENTREZID)])))}
entrez<-lapply(ens,function(s)lapply(s,map_entrez))
# Validate inputs against Venn workbooks before enrichment.
wb<-c(Up=file.path(repo,'results/historical','12dpd_all_res_sig_genes_by_stage_named_up_regulated.xlsx'),Down=file.path(repo,'results/historical','12dpd_all_res_sig_genes_by_stage_named_down_regulated.xlsx'))
validation<-list()
for(dir in names(loci))for(code in names(labels)){cond<-labels[[code]];h<-read_excel(wb[[dir]],sheet=paste('all',code),.name_repair='minimal');old<-sort(unique(as.character(na.omit(h$EnsemblID_Hs))));old<-old[nzchar(old)];new<-ens[[dir]][[cond]];validation[[paste(dir,code)]]<-data.frame(direction=dir,condition=cond,historical=length(old),reproduced=length(new),historical_only=length(setdiff(old,new)),reproduced_only=length(setdiff(new,old)),unique_ensembl=!anyDuplicated(new),unique_entrez=!anyDuplicated(entrez[[dir]][[cond]]))}
validation<-bind_rows(validation);write.csv(validation,file.path(out,'input_validation.csv'),row.names=FALSE);print(validation,row.names=FALSE);if(any(validation$historical_only|validation$reproduced_only))stop('Input mismatch')

# Compare three candidate WS40 TF pipelines against historical Tables S10/S12:
# A existing symbol-based reconstruction, B Ensembl target mapping, C Ensembl + historical localized batch correction.
suppressPackageStartupMessages({library(DESeq2);library(SummarizedExperiment);library(decoupleR);library(org.Hs.eg.db);library(limma);library(readxl)})
repo <- normalizePath(Sys.getenv('WS40_REPO', unset='.')); root <- file.path(repo,'data/reference')
out <- file.path(repo,'results/reproduced/collectri');dir.create(out,recursive=TRUE,showWarnings=FALSE)
se0 <- readRDS(file.path(repo,'data/input/WS40_12dpd_original_SE.rds'))
se0 <- se0[rowSums(assay(se0))>=5,,drop=FALSE]
short <- as.character(se0$treatment); full <- unname(c(t0='Normoxia',SA='4hA',LA='24hA',SR='2hR',LR='24hR')[short]);stopifnot(!anyNA(full))
ids <- as.character(SummarizedExperiment::colData(se0)$Sample);if(anyNA(ids)||anyDuplicated(ids))ids<-paste0(full,'_',ave(seq_along(full),full,FUN=seq_along))
colnames(se0)<-ids
key <- read.csv(file.path(repo,'data/input/Orthofinder_blast_merge_Hs_230628_MapIDS.csv'),fileEncoding='latin1',stringsAsFactors=FALSE,check.names=FALSE)
map <- data.frame(ensembl=sub('\\..*$','',as.character(key[[2]])),locus=as.character(key[[8]]),stringsAsFactors=FALSE)
map <- map[!is.na(map$ensembl)&nzchar(map$ensembl)&map$ensembl!='NULL'&!is.na(map$locus)&nzchar(map$locus),]
rawnet <- read.csv(file.path(root,'CollecTRI.csv'),stringsAsFactors=FALSE,check.names=FALSE)
net <- rawnet[,1:3];names(net)<-c('source','target','weight');net$weight<-as.numeric(net$weight);net<-net[!is.na(net$source)&nzchar(net$source)&!is.na(net$target)&nzchar(net$target)&is.finite(net$weight),]
net$target_ensembl <- unname(mapIds(org.Hs.eg.db,keys=net$target,column='ENSEMBL',keytype='SYMBOL',multiVals='first'))
net_ens <- net[!is.na(net$target_ensembl)&nzchar(net$target_ensembl),]
hist10 <- as.data.frame(readxl::read_excel(file.path(repo,'results/historical/Table_S10_Clouser_et_al_WS40_2026.xlsx'),.name_repair='minimal'))
hist12 <- as.data.frame(readxl::read_excel(file.path(repo,'results/historical/Table_S12_Clouser_et_al_WS40_2026.xlsx'),.name_repair='minimal'))
make_ens <- function(m){
 z<-m[rownames(m)%in%map$locus,,drop=FALSE];rownames(z)<-map$ensembl[match(rownames(z),map$locus)]
 a<-aggregate(z,by=list(rownames(z)),FUN=mean);rownames(a)<-a$Group.1;a$Group.1<-NULL;as.matrix(a)
}
score_variant <- function(label,mat){
 mat<-make_ens(mat)
 tf<-as.data.frame(run_wmean(mat=mat,network=net_ens,.source='source',.target='target_ensembl',.mor='weight',times=1000,seed=42))
 tf<-tf[tf$statistic=='norm_wmean',];tf$treatment<-full[match(tf$condition,ids)];stopifnot(!anyNA(tf$treatment))
 means<-aggregate(score~treatment+source,data=tf,mean);write.csv(means,file.path(out,paste0(label,'_condition_means.csv')),row.names=FALSE)
 fs<-do.call(rbind,lapply(split(tf,tf$source),function(z){q<-summary(aov(score~factor(treatment),data=z))[[1]];data.frame(source=z$source[1],F_statistic=q[['F value']][1],p_value=q[['Pr(>F)']][1])}));fs<-fs[order(-fs$F_statistic),];fs$rank<-seq_len(nrow(fs));write.csv(fs,file.path(out,paste0(label,'_global_F.csv')),row.names=FALSE)
 h10<-hist10[,c('treatment','source','mean_score')];h10$mean_score<-as.numeric(gsub(',', '', as.character(h10$mean_score)));j<-merge(h10,means,by=c('treatment','source'),suffixes=c('_historical','_reproduced'))
 h12<-hist12[,c('source','F_statistic')];h12$F_statistic<-as.numeric(h12$F_statistic);k<-merge(h12,fs,by='source',suffixes=c('_historical','_reproduced'))
 data.frame(variant=label,TFs=length(unique(tf$source)),S10_shared=length(unique(j$source)),S10_Pearson=cor(as.numeric(j$mean_score),as.numeric(j$score)),S10_Spearman=cor(as.numeric(j$mean_score),as.numeric(j$score),method='spearman'),S10_RMSE=sqrt(mean((as.numeric(j$score)-as.numeric(j$mean_score))^2)),S12_shared=nrow(k),S12_Pearson=cor(k$F_statistic_historical,k$F_statistic_reproduced),S12_Spearman=cor(k$F_statistic_historical,k$F_statistic_reproduced,method='spearman'),S12_RMSE=sqrt(mean((k$F_statistic_reproduced-k$F_statistic_historical)^2)),HSF1_F=fs$F_statistic[fs$source=='HSF1'],HSF2_F=fs$F_statistic[fs$source=='HSF2'],ERF_F=fs$F_statistic[fs$source=='ERF'])
}
# B: Ensembl, no batch correction.
se <- se0;se$treatment<-factor(short,levels=c('t0','SA','LA','SR','LR'));dds<-DESeqDataSet(se,design=~treatment);m1<-assay(vst(dds,blind=FALSE));m1<-sweep(m1,1,rowMeans(m1),'-')
ans <- score_variant('Ensembl_no_batch',m1)
write.csv(ans,file.path(out,'validated_no_batch_summary.csv'),row.names=FALSE)
print(ans,row.names=FALSE)

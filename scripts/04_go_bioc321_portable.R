suppressPackageStartupMessages(library(BiocManager))
needed <- c('clusterProfiler','org.Hs.eg.db','readxl','GOSemSim')
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) BiocManager::install(missing,ask=FALSE,update=FALSE)
suppressPackageStartupMessages({library(clusterProfiler);library(org.Hs.eg.db);library(readxl)})
out <- '/work/repo/results/reproduced/go_bioc321';dir.create(out,recursive=TRUE,showWarnings=FALSE)
base <- '/work/repo/results/historical'
refs <- '/work/repo/results/historical'
codes <- c(SA='4hA',LA='24hA',SR='2hR',LR='24hR')
getlists <- function(dir){p<-file.path(base,paste0('12dpd_all_res_sig_genes_by_stage_named_',tolower(dir),'_regulated.xlsx'));z<-lapply(names(codes),function(k){x<-readxl::read_excel(p,sheet=paste('all',k),.name_repair='minimal');v<-substr(as.character(x$EnsemblID_Hs),1,15);unique(v[!is.na(v)&nzchar(v)&v!='NULL'])});names(z)<-unname(codes);z[lengths(z)>0]}
key<-function(x)unique(paste(as.character(x$Cluster),as.character(x$ID),sep='|'))
rows<-list()
for(dir in c('Up','Down')){
 genes<-getlists(dir)
 h<-as.data.frame(readxl::read_excel(file.path(refs,if(dir=='Up')'Table_S7_Clouser_et_al_WS40_2026.xlsx' else 'Table_S9_Clouser_et_al_WS40_2026.xlsx')))
 for(mode in c('default','3_800')){
  a<-list(geneCluster=genes,fun='enrichGO',OrgDb=org.Hs.eg.db,keyType='ENSEMBL',ont='ALL',readable=TRUE)
  if(mode=='3_800')a<-c(a,list(minGSSize=3,maxGSSize=800,pvalueCutoff=.05,qvalueCutoff=.05))
  result<-do.call(clusterProfiler::compareCluster,a);d<-as.data.frame(result)
  write.csv(d,file.path(out,paste0(dir,'_',mode,'_raw.csv')),row.names=FALSE)
  x<-key(h);y<-key(d)
  rows[[paste(dir,mode)]]<-data.frame(release='3_21',direction=dir,setting=mode,historical=length(x),reproduced=length(y),shared=length(intersect(x,y)),historical_only=length(setdiff(x,y)),reproduced_only=length(setdiff(y,x)))
  if(dir=='Down')write.csv(d[d$Cluster=='24hA' & d$ID=='GO:1990204',],file.path(out,paste0('GO1990204_',mode,'.csv')),row.names=FALSE)
 }
}
s<-do.call(rbind,rows);write.csv(s,file.path(out,'comparison_summary.csv'),row.names=FALSE);print(s,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))

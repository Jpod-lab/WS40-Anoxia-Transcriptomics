#!/usr/bin/env Rscript
# Usage: Rscript scripts/validate_go_tables.R historical_S7.xlsx reproduced_Up_GO.csv historical_S9.xlsx reproduced_Down_GO.csv
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=4)stop('Expected historical S7, reproduced Up GO, historical S9, reproduced Down GO')
if(!requireNamespace('readxl',quietly=TRUE))stop('Install readxl')
for(i in c(1,3)){
 h<-as.data.frame(readxl::read_excel(args[i]));r<-read.csv(args[i+1],check.names=FALSE)
 k<-c('Cluster','ID');stopifnot(all(k%in%names(h)),all(k%in%names(r)))
 key<-function(x)paste(x$Cluster,x$ID,sep='|')
 a<-unique(key(h));b<-unique(key(r));cat('\n',basename(args[i]),': historical',length(a),'reproduced',length(b),'shared',length(intersect(a,b)),'\n')
 m<-merge(h,r,by=k,suffixes=c('_historical','_reproduced'))
 for(col in c('GeneRatio','BgRatio','Count','pvalue','p.adjust')){
  old<-paste0(col,'_historical');new<-paste0(col,'_reproduced')
  if(!all(c(old,new)%in%names(m)))next
  if(col%in%c('pvalue','p.adjust','Count')){
   x<-as.numeric(m[[old]]);y<-as.numeric(m[[new]]);cat(col,'max_abs_difference',max(abs(x-y),na.rm=TRUE),'\n')
  }else cat(col,'exact_rows',sum(as.character(m[[old]])==as.character(m[[new]]),na.rm=TRUE),'of',nrow(m),'\n')
 }
}

#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=4)stop('Usage: Rscript validate_go_complete.R S7.xlsx Up.csv S9.xlsx Down.csv')
suppressPackageStartupMessages(library(readxl))
for(i in c(1,3)){
 h<-as.data.frame(read_excel(args[i]));r<-read.csv(args[i+1],check.names=FALSE)
 m<-merge(h,r,by=c('Cluster','ID'),suffixes=c('_historical','_reproduced'))
 cat('\n',basename(args[i]),' historical=',nrow(h),' reproduced=',nrow(r),' matched=',nrow(m),'\n',sep='')
 for(field in c('GeneRatio','BgRatio','Count','geneID','pvalue','p.adjust')){
  a<-m[[paste0(field,'_historical')]];b<-m[[paste0(field,'_reproduced')]]
  if(is.null(a)||is.null(b))next
  if(field %in% c('pvalue','p.adjust')) cat(field,' max absolute difference: ',max(abs(as.numeric(a)-as.numeric(b)),na.rm=TRUE),'\n',sep='')
  else if(field=='geneID'){
   norm<-function(x) vapply(strsplit(as.character(x),'/',fixed=TRUE),function(v)paste(sort(unique(v)),collapse='/'),character(1))
   cat('geneID unordered-set matches: ',sum(norm(a)==norm(b)),'/',nrow(m),'\n',sep='')
  }else cat(field,' exact matches: ',sum(as.character(a)==as.character(b),na.rm=TRUE),'/',nrow(m),'\n',sep='')
 }
}

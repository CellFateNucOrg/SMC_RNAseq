library(ggplot2)
library(dplyr)
library(tidyr)
library(eulerr)
library(lattice)
library(gridExtra)
library(gplots)

source("functions.R")
source("./variableSettings.R")
scriptName <- "compareOrientation"
print(scriptName)
# note this script should be only used with the full genome's gene set

if(filterData){
  fileNamePrefix<-filterPrefix
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
} else {
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
}

makeDirs(outPath,dirNameList=paste0(c("plots/","txt/"),
                                    paste0(dirname(fileNamePrefix),"/",
                                           scriptName)))


# for gene-wise pileups
df<-NULL
for (grp in useContrasts[c(3,6,7,8,9)]){
  salmon<-GRanges(readRDS(paste0(outPath,"/rds/",fileNamePrefix,
                                 contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
  names(salmon)<-NULL
  if(is.null(df)){
    salmon1<-salmon
    mcols(salmon1)<-mcols(salmon1)[,c("wormbaseID","sequenceID","publicID","baseMean")]
    df<-data.frame(salmon1)
  }
  mcols(salmon)<-mcols(salmon)[,c("wormbaseID","log2FoldChange")]
  tmp<-data.frame(salmon)
  colnames(tmp)<-gsub("log2FoldChange",paste0("lfc_",grp),colnames(tmp))
  df<-left_join(df,tmp,by=c("seqnames","start","end","width","strand","wormbaseID" ))
}

sort(GRanges(df))
saveRDS(sort(GRanges(df)),file=paste0(outPath,"/rds/",fileNamePrefix,
                           "lfc_HiCsamples.rds"))



grList<-list()
grp=useContrasts[3]
for (grp in useContrasts[c(3,6,7,8,9)]){
  salmon<-GRanges(readRDS(paste0(outPath,"/rds/",fileNamePrefix,
                            contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
  names(salmon)<-NULL
  mcols(salmon)<-mcols(salmon)[,c("baseMean","log2FoldChange","padj","wormbaseID","publicID")]
  salmon$SMC<-grp
  salmon<-sort(salmon,ignore.strand=T)
  neighbours<-list()
  for(i in 1:length(salmon)){
    neighbours[i]<-nearestKNeighbors(salmon[i],salmon,ignore.strand=T,k=2)[[1]]
  }
}

ns<-salmon
strand(ns)<-"*"

gaps(salmon)
summary(width(gaps(ns)))

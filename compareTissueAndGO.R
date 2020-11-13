library(readxl)
library(ggVennDiagram)
library(ggplot2)
library(EnhancedVolcano)
library(wormcat)

source("functions.R")

outPath="."
padjVal=0.05
lfcVal=0.5

fileList<-read.table(paste0(outPath,"/fastqList.txt"),stringsAsFactors=F,header=T)

# extract the strain variable
strain<-factor(as.character(unique(fileList$sampleName),levels=c("366","382","775","784")))
SMC<-strain
levels(SMC)<-c("wt","dpy26cs","kle2cs","scc1cs")

controlGrp<-levels(SMC)[1] # control group
groupsOI<-levels(SMC)[-1]


############################
## Gene name lists for wormcat
############################

## all genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/salmon_",grp,"_DESeq2_fullResults.rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                          namePadjCol="padj",
                          nameLfcCol="log2FoldChange",
                          direction="both",
                          chr="all", nameChrCol="chr"))
}

sigGenes<-lapply(sigTables, "[", ,"wormbaseID")

for(i in 1:length(sigGenes)){
  write.table(sigGenes[[i]], file=paste0(outPath,"/txt/wormCat_",
                                         names(sigGenes)[i],"_padj",
                                         formatC(padjVal,format="e",digits=0),
                                         "_lfc", lfcVal,".txt"),
                                         quote=F,row.names=F, col.names=F)
}


### upregulated genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/salmon_",grp,"_DESeq2_fullResults.rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="gt",
                        chr="all", nameChrCol="chr"))
}


sigGenes<-lapply(sigTables, "[", ,"wormbaseID")

for(i in 1:length(sigGenes)){
  write.table(sigGenes[[i]], file=paste0(outPath,"/txt/wormCat_UP_",
                                         names(sigGenes)[i],"_padj",
                                         formatC(padjVal,format="e",digits=0),
                                         "_lfc", lfcVal,".txt"),
              quote=F,row.names=F, col.names=F)
}



### down regulated genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/salmon_",grp,"_DESeq2_fullResults.rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=-lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="lt",
                        chr="all", nameChrCol="chr"))
}


sigGenes<-lapply(sigTables, "[", ,"wormbaseID")

for(i in 1:length(sigGenes)){
  write.table(sigGenes[[i]], file=paste0(outPath,"/txt/wormCat_DOWN_",
                                         names(sigGenes)[i],"_padj",
                                         formatC(padjVal,format="e",digits=0),
                                         "_lfc", lfcVal,".txt"),
              quote=F,row.names=F, col.names=F)
}







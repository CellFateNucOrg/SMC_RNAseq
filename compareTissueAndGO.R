library(readxl)
library(ggVennDiagram)
library(ggplot2)
library(EnhancedVolcano)
library(wormcat)
library(xlsx)

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

  sigTables[[paste0(grp,"_all")]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                          namePadjCol="padj",
                          nameLfcCol="log2FoldChange",
                          direction="both",
                          chr="all", nameChrCol="chr"))
}
sigGenes<-lapply(sigTables, "[", ,"wormbaseID")




### upregulated genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/salmon_",grp,"_DESeq2_fullResults.rds"))

  sigTables[[paste0(grp,"_up")]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="gt",
                        chr="all", nameChrCol="chr"))
}
sigGenesUp<-lapply(sigTables, "[", ,"wormbaseID")



### down regulated genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/salmon_",grp,"_DESeq2_fullResults.rds"))

  sigTables[[paste0(grp,"_down")]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=-lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="lt",
                        chr="all", nameChrCol="chr"))

}

sigGenesDown<-lapply(sigTables, "[", ,"wormbaseID")


wormcatIn<-c(sigGenes,sigGenesUp,sigGenesDown)

wormcatIn<-lapply(wormcatIn,function (x) c("Wormbase ID",x))

openxlsx::write.xlsx(wormcatIn,
                     file=paste0(outPath,"/wormcat/wormcat.xlsx"))



#worm_cat_fun( file_to_process=paste0(outPath,"/wormcat/wormcat.xlsx"), output_dir=paste0(outPath,"/wormcat/"), input_type="Wormbase.ID")


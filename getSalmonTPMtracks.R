library(DESeq2)
library(AnnotationDbi)
library(Organism.dplyr)
library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
#library("TxDb.Celegans.UCSC.ce11.refGene")
#library("TxDb.Celegans.UCSC.ce11.ensGene")
library(tximport)
library(GenomicFeatures)
library(ggplot2)
library("RColorBrewer")
library("PoiClaClu")
library("pheatmap")
library(tidyr)
library(EnhancedVolcano)
library(affy)
library("gplots")
library(ggpubr)

#library(REBayes)
# get funciton for converting gene names from WBID to publicID
#source("~/Documents/MeisterLab/GenomeVer/geneNameConversion/convertingGeneNamesFunction1.R")
source("./variableSettings.R")
source("./functions.R")

scriptName <- "getSalmonTPM"
print(scriptName)

if(filterData){
  fileNamePrefix<-filterPrefix
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
} else {
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
}

makeDirs(outPath,dirNameList=paste0(c("plots/","tracks/","rds/"),
                                    paste0(dirname(fileNamePrefix),"/",
                                           scriptName)))



####
### some variables
#####

genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6], IRanges(start=1, end=seqlengths(Celegans)[1:6]))

wbseqinfo<-seqinfo(Celegans)
seqnames(wbseqinfo)<-c(gsub("chr","",seqnames(Celegans)))
seqnames(wbseqinfo)<-c(gsub("^M$","MtDNA",seqnames(wbseqinfo)))
genome(wbseqinfo)<-genomeVer
ce11seqinfo<-seqinfo(Celegans)



# Create metadata object --------------------------------------------------
###############################################################-
### create metadata ------
###############################################################-

if(!file.exists(paste0(outPath,"/wbGeneGR_WS275.rds"))){
  source("./createMetadataObj.R")
}

metadata<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))
#export.bed(metadata,paste0(outPath,"/wbGeneGR_WS275_ce11.bed"))

# load a txdb of wormbase data and create a tx2gene object
txdb<-loadDb(paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.", genomeVer,
                    ".annotations.sqlite"))
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")



###############################################################-
# Import into DESeq2 ------------------------------------------------------
###############################################################-

# salmon output: Name, length, estimated length (using sequence composition etc),
# TPM (read in as abundance in txi object) and numReads (read in as counts in txi object)
#
# when using different settings for "countsFromAbundance" it is the counts that
# change
#
# We don't really want double normalisation so use the TPM read in from salmon
# that are under teh abundance variable



# import the count matrices
txi<-tximport(sampleTable$fileName,type="salmon",tx2gene=tx2gene,
              countsFromAbundance="no")

 # txi1<-tximport(sampleTable$fileName,type="salmon",tx2gene=tx2gene,
 #               countsFromAbundance="scaledTPM")
 txi2<-tximport(sampleTable$fileName,type="salmon",tx2gene=tx2gene,
                 countsFromAbundance="lengthScaledTPM")

makeAvrTPMbedgraph<-function(strainName, txi, idx, metadata){
  gene<-data.frame("TPM"=rowMeans(txi$abundance[,idx]))
  gene$wormbaseID<-rownames(gene)

  md<-data.frame(metadata)
  gene<-left_join(gene,md,by="wormbaseID")

  genegr<-GRanges(seqnames=gene$seqnames,IRanges(gene$start,gene$end),gene$strand)
  genegr$score<-gene$TPM
  genegr$name<-gene$wormbaseID
  seqinfo(genegr)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)
  genegr<-sort(genegr)
  export.bed(genegr,con=paste0(outPath,"/tracks/PMW",strainName,"_TPM_avr.bed"))
  export(genegr,con=paste0(outPath,"/tracks/PMW",strainName,"_TPM_avr.bedgraph"),format="bedGraph")

  bins <- tileGenome(seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans), tilewidth=10,
                     cut.last.tile.in.chrom=TRUE)
  cov<-coverage(genegr,weight="score")
  ba<-binnedAverage(bins,numvar=cov,varname="score")
  export(ba,con=paste0(outPath,"/tracks/PMW",strainName,"_TPM_avr.bw"),format="bigwig")
}


getAvrTPM<-function(strainName, txi, idx, metadata){
  gene<-data.frame("TPM"=rowMeans(txi$abundance[,idx]))
  gene$wormbaseID<-rownames(gene)

  md<-data.frame(metadata)
  gene<-right_join(md,gene,by="wormbaseID")
  colnames(gene)<-gsub("TPM",paste0(strainName,"_TPM"),colnames(gene))
  return(gene)
}



options(tibble.width = Inf)
options(tibble.print_max = Inf)

tpmStrain="366"
idx<-sampleTable$strain==tpmStrain & sampleTable$auxin=="0mM"
sampleTable[idx,]  # for visual inspection
makeAvrTPMbedgraph(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)
tpm<-getAvrTPM(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)

tpmStrain="382"
idx<-sampleTable$strain==tpmStrain & sampleTable$auxin=="0mM"
sampleTable[idx,]  # for visual inspection
makeAvrTPMbedgraph(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)
x<-getAvrTPM(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)
tpm<-left_join(tpm, x)

tpmStrain="822"
idx<-sampleTable$strain==tpmStrain & sampleTable$auxin=="1mM"
sampleTable[idx,]  # for visual inspection
makeAvrTPMbedgraph(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)
x<-getAvrTPM(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)
tpm<-left_join(tpm, x)

tpmStrain="775"
idx<-sampleTable$strain==tpmStrain & sampleTable$auxin=="0mM"
sampleTable[idx,]  # for visual inspection
makeAvrTPMbedgraph(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)
x<-getAvrTPM(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)
tpm<-left_join(tpm, x)


tpmStrain="784"
idx<-sampleTable$strain==tpmStrain & sampleTable$auxin=="0mM"
sampleTable[idx,]  # for visual inspection
makeAvrTPMbedgraph(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)
x<-getAvrTPM(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)
tpm<-left_join(tpm, x)

tpmStrain="828"
idx<-sampleTable$strain==tpmStrain & sampleTable$auxin=="0mM"
sampleTable[idx,]  # for visual inspection
makeAvrTPMbedgraph(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)
x<-getAvrTPM(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)
tpm<-left_join(tpm, x)


tpmStrain="844"
idx<-sampleTable$strain==tpmStrain & sampleTable$auxin=="0mM"
sampleTable[idx,]  # for visual inspection
makeAvrTPMbedgraph(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)
x<-getAvrTPM(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)
tpm<-left_join(tpm, x)

tpmStrain="821"
idx<-sampleTable$strain==tpmStrain & sampleTable$auxin=="1mM"
sampleTable[idx,]  # for visual inspection
makeAvrTPMbedgraph(strainName=tpmStrain, txi=txi, idx=idx, metadata=metadata)

sort(GRanges(tpm))
saveRDS(sort(GRanges(tpm)),file=paste0(outPath,"/rds/",outputNamePrefix,
                                      "tpm_HiCsamples.rds"))




##########3
# get list of expressed genes
grp=useContrasts[1]
df<-NULL
for (grp in useContrasts[c(1,3,6,7,8,9)]){
  salmon<-data.frame(readRDS(paste0(outPath,"/rds/",fileNamePrefix,
                                 contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
  rownames(salmon)<-NULL
  idx<-!is.na(salmon$padj)
  salmon<-salmon[idx,c("chr","start","end", "strand",  "wormbaseID", "publicID", "sequenceID","baseMean","log2FoldChange")]
  colnames(salmon)<-gsub("log2FoldChange",paste0("lfc_",grp),colnames(salmon))
  if(is.null(df)){
    df<-salmon
  } else {
    df<-inner_join(df,salmon)
  }
  print(dim(df))
}

# add tpm
tpm<-data.frame(readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/2021_RNAseq_MDas/rds/p0.05_lfc0.5_filtChrAX/compareToChipSeq/filtChrAX_tpm_HiCsamples.rds"))
colnames(tpm)<-gsub("seqnames","chr",colnames(tpm))
df1<-left_join(df,tpm)
colnames(df1)<-gsub("X([[:digit:]]{3})_TPM","tpm_\\1",colnames(df1))
df1$width<-NULL
df1$entrezID<-NULL

gr<-sort(GRanges(df1))

saveRDS(gr,file=paste0(outPath,"/rds/",outputNamePrefix,
                             "lfc_tpm_expressedGenes_HiCsamples.rds"))

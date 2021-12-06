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

####
### some variables
#####

genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6], IRanges(start=1, end=seqlengths(Celegans)[1:6]))

wbseqinfo<-seqinfo(Celegans)
seqnames(wbseqinfo)<-c(gsub("chr","",seqnames(Celegans)))
seqnames(wbseqinfo)<-c(gsub("^M$","MtDNA",seqnames(wbseqinfo)))
genome(wbseqinfo)<-genomeVer
ce11seqinfo<-seqinfo(Celegans)

makeDirs(outPath,dirNameList=paste0(c("rds/","plots/","txt/","tracks/"),dirname(fileNamePrefix)))



# Create metadata object --------------------------------------------------
###############################################################-
### create metadata ------
###############################################################-

if(!file.exists(paste0(outPath,"/wbGeneGR_WS275.rds"))){
  source("./createMetadataObj.R")
}

metadata<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))


# load a txdb of wormbase data and create a tx2gene object
txdb<-loadDb(paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.", genomeVer,
                    ".annotations.sqlite"))
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")



###############################################################-
# Import into DESeq2 ------------------------------------------------------
###############################################################-


# import the count matrices
txi<-tximport(sampleTable$fileName,type="salmon",tx2gene=tx2gene,
              countsFromAbundance="no")

 # txi1<-tximport(sampleTable$fileName,type="salmon",tx2gene=tx2gene,
 #               countsFromAbundance="scaledTPM")
 # txi2<-tximport(sampleTable$fileName,type="salmon",tx2gene=tx2gene,
 #                countsFromAbundance="lengthScaledTPM")

idx366<-sampleTable$strain=="366"
gene<-data.frame("TPM"=rowMeans(txi$abundance[,idx366]))
gene$wormbaseID<-rownames(gene)

md<-data.frame(metadata)
gene<-left_join(gene,md,by="wormbaseID")

genegr<-GRanges(seqnames=gene$seqnames,IRanges(gene$start,gene$end),gene$strand)
genegr$score<-gene$TPM
genegr$name<-gene$wormbaseID
seqinfo(genegr)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)
genegr<-sort(genegr)
export.bed(genegr,con=paste0(outPath,"/tracks/PMW366_TPM_avr.bed"))
export(genegr,con=paste0(outPath,"/tracks/PMW366_TPM_avr.bedgraph"),format="bedGraph")


idx382<-sampleTable$strain=="382"
gene<-data.frame("TPM"=rowMeans(txi$abundance[,idx382]))
gene$wormbaseID<-rownames(gene)

md<-data.frame(metadata)
gene<-left_join(gene,md,by="wormbaseID")

genegr<-GRanges(seqnames=gene$seqnames,IRanges(gene$start,gene$end),gene$strand)
genegr$score<-gene$TPM
genegr$name<-gene$wormbaseID
seqinfo(genegr)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)
genegr<-sort(genegr)
export.bed(genegr,con=paste0(outPath,"/tracks/PMW382_TPM_avr.bed"))
export(genegr,con=paste0(outPath,"/tracks/PMW382_TPM_avr.bedgraph"),format="bedGraph")


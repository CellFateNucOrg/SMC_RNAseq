library(DESeq2)
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

source("~/Documents/MeisterLab/GenomeVer/geneNameConversion/convertingGeneNamesFunction1.R")
source("functions.R")
###############################################################
### some variables
###############################################################
fileNamePrefix="salmon_"
outPath="."
genomeVer="WS275"
genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)

genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6], IRanges(start=1, end=seqlengths(Celegans)[1:6]))

wbseqinfo<-seqinfo(Celegans)
seqnames(wbseqinfo)<-c(gsub("chr","",seqnames(Celegans)))
seqnames(wbseqinfo)<-c(gsub("^M$","MtDNA",seqnames(wbseqinfo)))
genome(wbseqinfo)<-genomeVer
ce11seqinfo<-seqinfo(Celegans)

makeDirs(outPath,dirNameList=c("rds","plots","txt","tracks"))


fileList<-read.table(paste0(outPath,"/fastqList.txt"),stringsAsFactors=F,header=T)


sampleNames<-paste(fileList$sampleName, fileList$repeatNum, fileList$laneNum, sep="_")

#fileNames<-paste0(outPath,"/salmon/ncRNA/",sampleNames,"/quant.sf")
#fileNames<-paste0(outPath,"/salmon/pseudoRNA/",sampleNames,"/quant.sf")
fileNames<-paste0(outPath,"/salmon/tnRNA/",sampleNames,"/quant.sf")


sampleTable<-data.frame(fileName=fileNames,sampleName=sampleNames,stringsAsFactors=F)

# extract the technical replicate variable
sampleTable$replicate=fileList$repeatNum
sampleTable$lane=fileList$laneNum

# extract the strain variable
sampleTable$strain<-factor(as.character(fileList$sampleName),levels=c("366","382","775","784"))
sampleTable$SMC<-sampleTable$strain
levels(sampleTable$SMC)<-c("wt","dpy26cs","kle2cs","scc1cs")

controlGrp<-levels(sampleTable$SMC)[1] # control group
groupsOI<-levels(sampleTable$SMC)[-1] # groups of interest to contrast to control

# metadata
metadata<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))

# load a txdb of wormbase data and create a tx2gene object
txdb<-loadDb(paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.", genomeVer,
                    ".annotations.sqlite"))
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")



###############################################################
### get samples into DESeq2
###############################################################

# import the count matrices
txi<-tximport(sampleTable$fileName,type="salmon",tx2gene=tx2gene)

# read samples into DESeq2
dds <- DESeqDataSetFromTximport(txi=txi,
                                colData=sampleTable,
                                design=~replicate+lane+SMC)



###############################################################
### DESeq2 differential expression analysis (using negative binomial distribution)
###############################################################

#dds<-collapseReplicates(dds,groupby=samples$sampleID,renameCols=T)
#dds <- DESeq(dds)
# This function performs a default analysis through the steps:
#   1. estimation of size factors: estimateSizeFactors
#   2. estimation of dispersion: estimateDispersions
#   3. Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
#   returns a DESeqDataSet object

idx<-match(rownames(dds),metadata$wormbaseID)
# add gene and chormosome names as metadata
featureData <- data.frame(gene=rownames(dds),
                          chr=as.vector(seqnames(metadata))[idx]) #,

rowData(dds) <- DataFrame(mcols(dds), featureData)

#only take expressed genes
#dds <- dds[ rowSums(counts(dds)) > 1, ]
#print("Number of expressed genes:")
#dim(dds)[1]
#16606 genes from 20127

dds<-DESeq(dds)

statsPerSample<-data.frame(t(apply(counts(dds),2,summary)))
statsPerSample$totalCounts<-colSums(counts(dds))
rownames(statsPerSample)<-colData(dds)$sampleName
colnames(statsPerSample)<-c("min", "Q1", "median", "mean", "Q3", "max","totalCounts")
statsPerSample$zeros <- apply(counts(dds)==0, 2, sum)
statsPerSample$percZeros <- round(100*statsPerSample$zeros/nrow(counts(dds)),1)
print(statsPerSample)

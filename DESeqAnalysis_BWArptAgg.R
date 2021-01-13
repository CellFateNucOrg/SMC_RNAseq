library(DESeq2)
library(Organism.dplyr)
library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
#library("TxDb.Celegans.UCSC.ce11.refGene")
#library("TxDb.Celegans.UCSC.ce11.ensGene")
library(tximport)
library(GenomicFeatures)
library(GenomeInfoDb)
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
source("~/Documents/MeisterLab/GenomeVer/geneNameConversion/convertingGeneNamesFunction1.R")
source("functions.R")
####
### some variables
####
plotPDFs=F
fileNamePrefix="BWA_"
filterPrefix="BWA_rptFam_random_"
dataset="_random_rptFam"
filterData=T


padjVal=0.05
lfcVal=0.5
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

fileList<-read.table(paste0(outPath,"/fastqList.txt"), stringsAsFactors=F, header=T)


sampleNames<-paste(fileList$sampleName, fileList$repeatNum, fileList$laneNum,
                   fileList$libType, sep="_")
fileNames<-paste0(outPath,"/htseq/",sampleNames,"_union",dataset,".txt")

sampleTable<-data.frame(fileName=fileNames,sampleName=sampleNames,stringsAsFactors=F)

sampleTable$replicate=fileList$repeatNum
sampleTable$lane=fileList$laneNum

# extract the strain variable
sampleTable$strain<-factor(as.character(fileList$sampleName),levels=c("366","382","775","784"))
sampleTable$SMC<-sampleTable$strain
levels(sampleTable$SMC)<-c("wt","dpy26cs","kle2cs","scc1cs")

controlGrp<-levels(sampleTable$SMC)[1] # control group
groupsOI<-levels(sampleTable$SMC)[-1] # groups of interest to contrast to control


# Create metadata object --------------------------------------------------
###############################################################-
### create metadata
###############################################################-

if(!file.exists(paste0(outPath,"/wbGeneGRandRpts_WS275.rds"))){
  source(paste0(outPath,"/createMetadaObj.R"))
}

#metadata<-readRDS(paste0(outPath,"/wbGeneGRandRpts_WS275.rds"))
metadata<-readRDS(paste0(outPath,"metadataTbl_genes-rpts.rds"))


###############################################################-
# Import into DESeq2 ------------------------------------------------------
###############################################################-

sampleTable1<-sampleTable
sampleTable1$sampleName<-basename(sampleTable$fileName)

# read samples into DESeq2
dds <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable1,
                    directory=dirname(sampleTable1$fileName[1]),
                                design=~lane+SMC)

colData(dds)$sampleName<-paste(gsub(paste0("_union", dataset,
                                ".txt"), "",
                                sampleTable1$sampleName),
                               sep="_")

###############################################################-
### DESeq2 differential expression analysis (using negative binomial distribution)
###############################################################-

#dds<-collapseReplicates(dds,groupby=samples$sampleID,renameCols=T)
#dds <- DESeq(dds)
# This function performs a default analysis through the steps:
#   1. estimation of size factors: estimateSizeFactors
#   2. estimation of dispersion: estimateDispersions
#   3. Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
#   returns a DESeqDataSet object

rownames(dds)<-gsub("Gene:","",rownames(dds))

idx<-match(rownames(dds),metadata$ID)
# add gene and chormosome names as metadata
featureData <- data.frame(gene=rownames(dds),
                          chr=as.vector(metadata$seqnames)[idx],
                          bioType=metadata$bioType[idx])

rowData(dds) <- DataFrame(mcols(dds), featureData)

dds<-DESeq(dds)

# remove non-repeat genes
if(filterData){
  dds<-dds[grep("repeatFam",rowData(dds)$bioType),]
  fileNamePrefix<-filterPrefix
}

######################################################-
# Basic sample stats ------------------------------------------------------
######################################################-

## basic sample stats
sink(file=paste0(outPath,"/txt/", fileNamePrefix,
                 "all_logfile.txt"),
     append=FALSE, type="output")
statsPerSample<-data.frame(t(apply(counts(dds),2,summary)))
statsPerSample$totalCounts<-colSums(counts(dds))
rownames(statsPerSample)<-colData(dds)$sampleName
colnames(statsPerSample)<-c("min", "Q1", "median", "mean", "Q3", "max","totalCounts")
statsPerSample$zeros <- apply(counts(dds)==0, 2, sum)
statsPerSample$percZeros <- round(100*statsPerSample$zeros/nrow(counts(dds)),1)
print(statsPerSample)
sink()

pdf(file=paste0(outPath,"/plots/",fileNamePrefix,
                "sampleQC.pdf"), width=8,height=8,paper="a4")

#######-
# sample counts summary: boxplots and density plots -----------------------
#######-
## box plots
epsilon <- 1 # pseudo-count to avoid problems with log(0)
#df<-as.data.frame(log2(counts(dds) + epsilon))
#colnames(df)<-colData(dds)$sampleName
#ldf<-gather(df,key=sampleName,value=log2Counts)
#ldf$strain<-gsub("_.*$","",ldf$sampleName)
#p<-ggplot(ldf,aes(x=logFC,y=sampleName))+geom_boxplot(aes(fill=strain))
par(mar=c(5,7,4,2))
boxplot(log2(counts(dds) + epsilon), col=c("grey","blue","red","darkgreen")[as.factor(colData(dds)$strain)], pch=".",
        horizontal=TRUE, cex.axis=1,
        las=1, ylab=NULL, names=colData(dds)$sampleName,
        xlab=paste0("log2(counts ","+",epsilon,")"))
par(mar=c(5,4,4,2))

## density plots
plotDensity(log2(counts(dds) + epsilon), lty=1, col=as.factor(colData(dds)$sampleName), lwd=2, xlab="log2 counts per gene")
grid()
legend("topright", legend=colData(dds)$sampleName, col=as.factor(colData(dds)$sampleName), lwd=2,cex=0.8)



#########-
# dispersion estimates ----------------------------------------------------
#########-
plotDispEsts(dds)


#########-
# sample to sample heatmap ------------------------------------------------
#########-
#vsd <- vst(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds,blind=TRUE)
colnames(vsd)<-colData(dds)$sampleName
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colData(dds)$sampleName
colnames(sampleDistMatrix) <- colData(dds)$sampleName
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#########-
# Heatmap - most highly expressed genes -----------------------------------
#########-
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:184]
df <- as.data.frame(colData(dds)[,c("strain","replicate","lane")])
rownames(df)<-colData(dds)$sampleName
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df, main="184 repeat famlies - no clustering")

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=TRUE, annotation_col=df, main="184 repeat famlies - clustered by columns")



###########-
# pca ---------------------------------------------------------------------
###########-
p1<-plotPCA(vsd, intgroup=c("strain"))
print(p1)
p2<-plotPCA(vsd, intgroup=c("replicate"))
print(p2)
p3<-plotPCA(vsd, intgroup=c("lane"))
print(p3)
p4<-plotPCA(vsd, intgroup=c("strain", "replicate"))
print(p4)
dev.off()



##########-
# pairwise correlation between genes in different replicates --------------
##########-
#Define a function to draw a scatter plot for a pair of variables (samples) with density colors
plotFun <- function(x,y){
  dns <- densCols(x,y);
  points(x,y, col=dns, pch=".", panel.first=grid());
  #  abline(a=0, b=1, col="brown")
}

corFun <- function(x,y){
  par(usr = c(0, 1, 0, 1))
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 1.5*cor(x, y))
}

df<-as.data.frame(log2(counts(dds) + epsilon))
colnames(df)<-colData(dds)$sampleName
for (grp in c(controlGrp,groupsOI)){
  if(plotPDFs==T){
    pdf(file=paste0(outPath,"/plots/",fileNamePrefix, grp, "_correlations.pdf"),
        width=8,height=8,units="in",paper="a4")
  } else {
    png(file=paste0(outPath,"/plots/",fileNamePrefix, grp, "_correlations.png"),
        width=8,height=8,units="in",res=150)
  }
  idx<-colData(dds)$SMC %in% c(grp)
  pairs(df[,idx], panel=plotFun, lower.panel=corFun, labels=colnames(df)[idx], main=grp)
  dev.off()
}


##############################################################-
# Significant genes -------------------------------------------------------
##############################################################-


#res=list()
#resLFC=list()


for(grp in groupsOI){
  res<-results(dds,contrast=c("SMC",grp,controlGrp))
  sink(file=paste0(outPath,"/txt/",fileNamePrefix, grp,
                   "_logfile.txt"), append=TRUE,
       type="output")
  cat(paste0("Number of genes that change expression in ",grp," at different padj cutoffs:\n"))
  print(summary(res))
  print(summary(res,alpha=0.05))
  print(summary(res,alpha=0.01))
  sink()

  ### add metadata
  res$ID<-rownames(res)

  # shrink LFC estimates
  #resultsNames(dds) # to get names of coefficients
  resLFC<-lfcShrink(dds,coef=paste0("SMC_",grp,"_vs_",controlGrp), type="apeglm", res=res)
  class(resLFC)
  ### add metadata
  resLFC$ID<-rownames(resLFC)
  idx<-match(rownames(resLFC),metadata$ID)

  resLFC$rptID<-as.vector(metadata$rptID)[idx]
  resLFC$rptfamName<-as.vector(metadata$rptfamName)[idx]
  resLFC$rptfamID<-as.vector(metadata$rptfamID)[idx]
  resLFC$rptType<-as.vector(metadata$rptType)[idx]
  resLFC$famSize<-as.vector(metadata$famSize)[idx]

  saveRDS(resLFC,file=paste0(outPath,"/rds/", fileNamePrefix, grp,
                             "_DESeq2_fullResults.rds"))

  #export csv with ordered results
  write.csv(resLFC[order(resLFC$padj),],
            file=paste0(outPath,"/txt/", fileNamePrefix,grp,
                        "_DESeq2_resultsTable.csv"),
            quote=F,row.names=F)

  # remove NAs
  res<-res[! is.na(res$padj),]
  resLFC<-resLFC[! is.na(resLFC$padj),]


  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,grp,
                  "_hclust_mostChanged.pdf"), width=8,height=11,paper="a4")

  #######-
  # plot results filtering threshold ----------------------------------------
  #######-
  plot(metadata(resLFC)$filterNumRej,
       type="b", ylab="number of rejections",
       xlab="quantiles of filter",main="Threshold for independant filtering of results")
  lines(metadata(resLFC)$lo.fit, col="red")
  abline(v=metadata(resLFC)$filterTheta)
  legend("topright",legend=paste0("Read count \nthreshold: ",
                                  round(metadata(resLFC)$filterThreshold,2)))


  ##########-
  # heirarchical clustering of most significantly changed genes -------------
  ##########-
  # select gene names based on FDR (5%)
  gene.kept <- rownames(resLFC)[resLFC$padj <= padjVal & !is.na(resLFC$padj) & abs(resLFC$log2FoldChange)>lfcVal]

  # Retrieve the normalized counts for gene of interest
  countTable.kept <- log2(counts(dds) + epsilon)[gene.kept, ]
  if(all(!is.null(dim(countTable.kept)), nrow(countTable.kept>2))){
    colnames(countTable.kept)<-colData(dds)$sampleName

    # Perform the hierarchical clustering with
    # A distance based on Pearson-correlation coefficient
    # and average linkage clustering as agglomeration criteria
    heatmap.2(as.matrix(countTable.kept),
            scale="row",
            hclust=function(x) hclust(x,method="average"),
            distfun=function(x) as.dist((1-cor(t(x)))/2),
            margin=c(6,0),
            trace="none",
            density="none",
            labRow="",
            #labCol = names(countTable.kept),
            cexCol=1,
            main=paste0(grp," changed genes (p<",padjVal,", lfc>0.5)"))
  }

  dev.off()


  ##########-
  # plot individual genes ---------------------------------------------------
  ##########-

  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,grp,
                  "_topGenes_normCounts.pdf"), width=8,height=11,paper="a4")
  par(mfrow=c(5,2))
  selectedGenes <- rownames(resLFC)[order(resLFC$padj)][1:20]
  for (g in selectedGenes) {
    barplot(counts(dds, normalized=TRUE)[g,],
            col=colData(dds)$SMC,
            main=resLFC$rptfamName[resLFC$ID==g], las=2, cex.names=0.8,
            names.arg=colData(dds)$sampleName)
    legend("topleft",legend=levels(colData(dds)$SMC),fill=c(1,2), cex=0.7)
  }
  dev.off()



  # MAplots -----------------------------------------------------------------
  ###########-
  # MAplot ALL genes
  ############-

  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,grp,
                  "_MAplots_results.pdf"), width=5,height=5,paper="a4")


  plotMA(res, main=paste0(grp," uncorrected LFC, threshold=", padjVal), ylim=c(-3,3), alpha=padjVal)
  plotMA(resLFC, main=paste0(grp," apeglm shrunk LFC, threshold=", padjVal), ylim=c(-3,3), alpha=padjVal)

  dev.off()




  #############-
  # Volcano plots -----------------------------------------------------------
  #############-
  #https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
  #all black plot for manual changing of colours.
  #pdf(file=paste0(outPath,"/plots/",fileNamePrefix, grp,
  #                "_volcanoPlot_allGenes.pdf"), width=8,height=6,paper="a4")

  keyvals<-rep('black', nrow(resLFC))
  names(keyvals)<-rep('NS',nrow(resLFC))
  keyvals[which(resLFC$padj<padjVal & abs(resLFC$log2FoldChange)> lfcVal)]<-'red'
  names(keyvals)[which(resLFC$padj<padjVal & abs(resLFC$log2FoldChange)> lfcVal)]<-paste0('p<',padjVal,' |lfc|>',lfcVal)
  sigUp<-sum(resLFC$padj<padjVal & resLFC$log2FoldChange> lfcVal)
  sigDown<-sum(resLFC$padj<padjVal & resLFC$log2FoldChange< -lfcVal)
  p1<-EnhancedVolcano(resLFC,
                      lab=paste0(resLFC$rptfamName,"(",
                                 resLFC$famSize, ")"),
                      labSize=3,
                      #labCol="#11111100",
                      labCol="black",
                      drawConnectors=T,
                      x="log2FoldChange",
                      y="padj",
                      #selectLab=rownames(resLFC)[12366],
                      #xlim=c(-5.5,5.5),
                      #ylim=c(0,80),
                      title= paste0(gsub("_$","",fileNamePrefix)," ",grp,
                                    " vs ", controlGrp),
                      titleLabSize=14,
                      subtitle=NULL,
                      caption = paste0(nrow(resLFC), ' expressed families. ',
                                       sigUp, " up, ",sigDown," down."),
                      captionLabSize = 12,
                      pCutoff=padjVal,
                      FCcutoff=lfcVal,
                      xlab=bquote(~Log[2]~'fold change'~.(grp)~'/'~.(controlGrp)),
                      ylab=bquote(~-Log[10]~adjusted~italic(P)),
                      #.legend=c('NS','P & Log2 FC'),
                      #legendLabels=c('NS', expression(p-value<padjVal~and~log[2]~FC>1)),
                      legendPosition = 'top',
                      legendLabSize = 12,
                      legendIconSize = 3.0,
                      axisLabSize=14,
                      colCustom=keyvals,
                      #col = c("black", "red"),
                      colAlpha=0.5,
                      pointSize = 1.0)
  #dev.off()
  if(plotPDFs==T){
    ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                           "_volcanoPlot_allGenes.pdf"), plot=p1,
           device="pdf",path=outPath, width=12,height=12,units="cm")
  } else {
    ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                           "_volcanoPlot_allGenes.png"), plot=p1,
           device="png",path=outPath, width=12,height=12,units="cm")
  }


  ##### create filtered tables of gene names
  results<-readRDS(paste0(outPath,"/rds/",fileNamePrefix, grp,
                          "_DESeq2_fullResults.rds"))
  #results<-na.omit(results)


  nrow(filterResults(results,padj=0.05,lfc=0,"both","all", writeTable=F))
  nrow(filterResults(results,padj=0.05,lfc=0.5,"both","all", writeTable=F))
  nrow(filterResults(results,padj=0.05,lfc=0.75,"both","all", writeTable=F))
  nrow(filterResults(results,padj=0.05,lfc=1,"both","all", writeTable=F))
  nrow(filterResults(results,padj=0.01,lfc=0,"both","all", writeTable=F))
  nrow(filterResults(results,padj=0.01,lfc=0.5,"both","all", writeTable=F))
  nrow(filterResults(results,padj=0.01,lfc=0.75,"both","all", writeTable=F))
  nrow(filterResults(results,padj=0.01,lfc=1,"both","all", writeTable=F))

}



# # Volcano - colour by other datasets --------------------------------------
#
#
#
# oscillating<-read.delim(paste0(outPath,"/oscillatingGenes.tsv"),header=T,
#                         stringsAsFactors=F)
#
# hsUp<-readRDS("hsUp_garrigues2019.rds")
# hsDown<-readRDS("hsDown_garrigues2019.rds")
#
#
# #grp=groupsOI[3]
#
# for(grp in groupsOI){
#   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))
#
#   #### oscillating genes
#   keyvals<-rep('black', nrow(salmon))
#   names(keyvals)<-rep('Other',nrow(salmon))
#   idx<-salmon$ID %in% oscillating$WB_ID
#   keyvals[idx]<-'red'
#   names(keyvals)[idx]<-"Oscillating"
#   sigUp<-sum(rowSums(cbind(salmon$padj< padjVal,
#                            salmon$log2FoldChange>lfcVal,
#                            keyvals=="red"), na.rm=T)==3, na.rm=T)
#   sigDown<-sum(rowSums(cbind(salmon$padj< padjVal,
#                              salmon$log2FoldChange<lfcVal,
#                              keyvals=="red"), na.rm=T)==3, na.rm=T)
#   p1<-EnhancedVolcano(salmon,
#                       lab=salmon$publicID,
#                       labSize=0.5,
#                       labCol="#11111100",
#                       x="log2FoldChange",
#                       y="padj",
#                       selectLab=salmon$publicID[12366],
#                       xlim=c(-5.5,5.5),
#                       ylim=c(0,65),
#                       title= paste0(grp," vs ", controlGrp),
#                       subtitle=NULL,
#                       caption = paste0(sum(keyvals=="red"), ' oscillating genes (Meeuse 2020). ',sigUp, " up, ",sigDown," down."),
#                       captionLabSize = 12,
#                       pCutoff=padjVal,
#                       FCcutoff=lfcVal,
#                       xlab=bquote(~Log[2]~'fold change'~.(grp)~'/'~.(controlGrp)),
#                       ylab=bquote(~-Log[10]~adjusted~italic(P)),
#                       #.legend=c('NS','P & Log2 FC'),
#                       #legendLabels=c('NS', expression(p-value<padjVal~and~log[2]~FC>1)),
#                       legendPosition = 'top',
#                       legendLabSize = 12,
#                       legendIconSize = 3.0,
#                       axisLabSize=14,
#                       colCustom=keyvals,
#                       #col = c("black", "red"),
#                       colAlpha=0.5,
#                       pointSize = 1.0)
#   #dev.off()
#   if(plotPDFs==T){
#     ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
#                            "_volcanoPlot_oscillatingGenes.pdf"), plot=p1,
#            device="pdf",path=outPath, width=12,height=12,units="cm")
#   } else {
#     ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
#                            "_volcanoPlot_oscillatingGenes.png"), plot=p1,
#            device="png",path=outPath, width=12,height=12,units="cm")
#   }
#
#
#
#   #### heat shock genes
#   myCols<-c("#11111100","red","red") # background, dataset1, dataset2
#   keyvals<-rep(myCols[1], nrow(salmon))
#   names(keyvals)<-rep('Other',nrow(salmon))
#   idx<-salmon$ID %in% hsUp$WormBase.ID | salmon$ID %in% hsDown$WormBase.ID
#   keyvals[idx]<-myCols[2]
#   names(keyvals)[idx]<-"Heatshock"
#   sigUp<-sum(rowSums(cbind(salmon$padj< padjVal,
#                            salmon$log2FoldChange>lfcVal,
#                            keyvals==myCols[2]), na.rm=T)==3, na.rm=T)
#   sigDown<-sum(rowSums(cbind(salmon$padj< padjVal,
#                              salmon$log2FoldChange< -lfcVal,
#                              keyvals==myCols[2]), na.rm=T)==3, na.rm=T)
#   p1<-EnhancedVolcano(salmon,
#                       lab=salmon$publicID,
#                       labSize=0.5,
#                       labCol=myCols[1],
#                       x="log2FoldChange",
#                       y="padj",
#                       selectLab=salmon$publicID[12366],
#                       xlim=c(-5.5,5.5),
#                       ylim=c(0,65),
#                       title= paste0(grp," vs ", controlGrp),
#                       subtitle=NULL,
#                       caption = paste0(sum(keyvals==myCols[2]), ' heatshock genes (Garrigues 2019). ',sigUp, " up, ",sigDown," down."),
#                       captionLabSize = 12,
#                       pCutoff=padjVal,
#                       FCcutoff=lfcVal,
#                       xlab=bquote(~Log[2]~'fold change'~.(grp)~'/'~.(controlGrp)),
#                       ylab=bquote(~-Log[10]~adjusted~italic(P)),
#                       #.legend=c('NS','P & Log2 FC'),
#                       #legendLabels=c('NS', expression(p-value<padjVal~and~log[2]~FC>1)),
#                       legendPosition = 'top',
#                       legendLabSize = 12,
#                       legendIconSize = 3.0,
#                       axisLabSize=14,
#                       colCustom=keyvals,
#                       #col = c("black", "red"),
#                       colAlpha=0.5,
#                       pointSize = 1.0)
#   #dev.off()
#   if(plotPDFs==T){
#     ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
#                            "_volcanoPlot_hsGenes.pdf"), plot=p1,
#            device="pdf",path=outPath, width=12,height=12,units="cm")
#   } else {
#     ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
#                            "_volcanoPlot_hsGenes.png"), plot=p1,
#            device="png",path=outPath, width=12,height=12,units="cm")
#   }
# }

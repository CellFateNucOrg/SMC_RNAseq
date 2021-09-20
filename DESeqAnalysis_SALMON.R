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

makeDirs(outPath,dirNameList=paste0(c("rds/","plots/","txt/","tracks/"),paste0("p",padjVal,"_lfc",lfcVal)))



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
txi<-tximport(sampleTable$fileName,type="salmon",tx2gene=tx2gene)

# read samples into DESeq2
dds <- DESeqDataSetFromTximport(txi=txi,
                                colData=sampleTable,
                          design = formula(modelTxt))



###############################################################-
### DESeq2 differential expression analysis (using negative binomial distribution)-----
###############################################################-

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

#remove unmapped or mitochondrial genes
idx<-is.na(rowData(dds)$chr) | rowData(dds)$chr %in% c("MtDNA","chrM")
dds<-dds[!idx,]

#prefilter rows with less than 10 reads in total
dds<-dds[rowSums(counts(dds)) >= 10,]
print(paste0(dim(dds)[1], " genes with >=10 reads total"))

###################-
####### filter genes-------------------------
#######-###########-
if(filterData){
   # remove filtered genes
   idx<-rowData(dds)$gene %in% toFilter
   dds<-dds[!idx,]

   fileNamePrefix=filterPrefix
}

#dds<-DESeq(dds)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit=1000)

# head(model.matrix(design(dds),colData(dds)))
# resultsNames(dds)
colMeans(coef(dds))

saveRDS(dds,file=paste0(outPath,"/rds/dds_object.rds"))

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
#explanation of this graph.
#https://github.com/hbctraining/DGE_workshop/blob/master/lessons/04_DGE_DESeq2_analysis.md

#########-
# sample to sample heatmap ------------------------------------------------
#########-
vsd <- vst(dds, blind=TRUE)
colnames(vsd)<-colData(dds)$sampleName
sampleDists <- stats::dist(t(assay(vsd)))
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
                decreasing=TRUE)[1:500]
df <- as.data.frame(colData(dds)[,c("strain","TIR1","auxin")])
rownames(df)<-colData(dds)$sampleName
#pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df, main="Top 500 expressed genes - no clustering")

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=TRUE, annotation_col=df, main="Top 500 expressed genes - clustered by columns")



###########-
# pca ---------------------------------------------------------------------
###########-
# get variables about data (ignore file and sample names)
colourByVar<-colnames(sampleTable)[grep("fileName|sampleName",colnames(sampleTable),invert=T)]
numPages<-length(colourByVar)/2
close
plotlist<-list()
for(varToUse in colourByVar){
   p1<-plotPCA(vsd, intgroup=c(varToUse),ntop=1000)+ggtitle(varToUse)
   plotlist[[varToUse]]<-p1
}
p<-ggarrange(plotlist=plotlist,ncol=1,nrow=3)
print(p)
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
          width=8,height=8,paper="a4")
   } else {
      png(file=paste0(outPath,"/plots/",fileNamePrefix, grp, "_correlations.png"),
       width=8,height=8,units="in",res=150)
   }
   idx<-colData(dds)[,varOI] %in% c(grp)
   pairs(df[,idx], panel=plotFun, lower.panel=corFun, labels=colnames(df)[idx], main=grp)
   dev.off()
}


##############################################################-
# Significant genes -------------------------------------------------------
##############################################################-
#res=list()
#resLFC=list()
for(grp in names(contrastNames)){
   ###############-
   # filter at the predetermined lfcVal -------
   # #############-
   if(is.character(contrastsOI[[contrastNames[[grp]]]])){
      res<-results(dds,contrast=list(contrastNames[[grp]]),alpha=padjVal)
   } else {
      res<-results(dds,contrast=contrastsOI[[contrastNames[[grp]]]],alpha=padjVal)
   }

   pdf(file=paste0(outPath,"/plots/",fileNamePrefix,contrastNames[[grp]],"_independentFilter_",
                   padjVal,".pdf"), width=8, height=8, paper="a4")
   plot(metadata(res)$filterNumRej,
        type="b", ylab="number of rejections",
        xlab="quantiles of filter",
        main=paste0("Independant filtering, ",grp,", alpha=",padjVal))
   lines(metadata(res)$lo.fit, col="red")
   abline(v=metadata(res)$filterTheta)
   legend("topright",legend=paste0("Mean norm count \nthreshold: ", round(metadata(res)$filterThreshold,2)))
   dev.off()
   # shrink LFC estimates
   #resultsNames(dds) # to get names of coefficients

   if(is.character(contrastsOI[[contrastNames[[grp]]]])){
      resLFC<-lfcShrink(dds,type="ashr",res=res)
   } else {
      resLFC<-lfcShrink(dds,type="ashr",res=res)
   }

   #resLFC<-lfcShrink(dds,coef=grp, type="apeglm", res=res)
   class(resLFC)
   ### add metadata
   resLFC$wormbaseID<-rownames(resLFC)
   idx<-match(rownames(resLFC),metadata$wormbaseID)
   resLFC$chr<-factor(seqnames(metadata),levels=paste0("chr",c("I","II","III","IV","V","X")))[idx]
   resLFC$start<-as.vector(start(metadata))[idx]
   resLFC$end<-as.vector(end(metadata))[idx]
   resLFC$strand<-as.vector(strand(metadata))[idx]
   resLFC$publicID<-as.vector(metadata$publicID)[idx]
   resLFC$sequenceID<-as.vector(metadata$sequenceID)[idx]
   resLFC$entrezID<-as.vector(metadata$entrezID)[idx]


   saveRDS(resLFC,file=paste0(outPath,"/rds/", fileNamePrefix, contrastNames[[grp]],
                              "_DESeq2_fullResults_p",padjVal,".rds"))

   #export csv with ordered results
   write.csv(resLFC[order(resLFC$padj),],
             file=paste0(outPath,"/txt/", fileNamePrefix,contrastNames[[grp]],
                         "_DESeq2_resultsTable_p",padjVal,".csv"),
             quote=F,row.names=F)

   # remove NAs from chr (unmapped or mtDNA) and padj (below filter threshold) columns
   idx<-is.na(resLFC$chr) | is.na(resLFC$padj)
   res<-res[!idx,]
   resLFC<-resLFC[!idx,]


   pdf(file=paste0(outPath,"/plots/",fileNamePrefix,contrastNames[[grp]],
                   "_hclust_mostChanged.pdf"), width=8,height=11,paper="a4")



   ##########-
   # heirarchical clustering of most significantly changed genes -------------
   ##########-
   # select gene names based on FDR (5%)
   gene.kept <- rownames(resLFC)[resLFC$padj <= padjVal & !is.na(resLFC$padj) & abs(resLFC$log2FoldChange)>lfcVal]

   # Retrieve the normalized counts for gene of interest
   countTable.kept <- log2(counts(dds) + epsilon)[gene.kept, ]
   dim(countTable.kept)
   colnames(countTable.kept)<-colData(dds)$sampleName

   # Perform the hierarchical clustering with
   # A distance based on Pearson-correlation coefficient
   # and average linkage clustering as agglomeration criteria
   heatmap.2(as.matrix(countTable.kept),
             scale="row",
             hclust=function(x) stats::hclust(x,method="average"),
             distfun=function(x) stats::as.dist((1-cor(t(x)))/2),
             margin=c(6,0),
             trace="none",
             density="none",
             labRow="",
             #labCol = names(countTable.kept),
             cexCol=1,
             main=paste0(grp," changed genes (p<",padjVal,", lfc>",lfcVal,")"))

   dev.off()


   ##########-
   # plot individual genes ---------------------------------------------------
   ##########-

   pdf(file=paste0(outPath,"/plots/",fileNamePrefix,contrastNames[[grp]],
                   "_topGenes_normCounts.pdf"), width=8,height=11,paper="a4")
   par(mfrow=c(5,2))
   selectedGenes <- rownames(resLFC)[order(resLFC$padj)][1:20]
   for (g in selectedGenes) {
      barplot(counts(dds, normalized=TRUE)[g,],
              col=colData(dds)[,varOI],
              main=g, las=2, cex.names=1,names.arg=colData(dds)$sampleName)
      legend("topleft",legend=levels(colData(dds)[,varOI]),fill=c(1,2), cex=0.7)
   }
   dev.off()


   ##########-
   # make GRanges for LFC ----------------------------------------------------
   ##########-
   #remove nas
   resGR<-GenomicRanges::GRanges(seqnames=resLFC$chr,
                                 IRanges::IRanges(start=resLFC$start,
                                                  end=resLFC$end),
                                 strand=resLFC$strand)
   seqlengths(resGR)<-seqlengths(Celegans)[1:6]
   mcols(resGR)<-resLFC[,c("wormbaseID","log2FoldChange","padj")]

   names(mcols(resGR))[names(mcols(resGR))=="log2FoldChange"]<-"score"
   resGR<-sort(resGR,ignore.strand=TRUE)

   #https://github.com/hochwagenlab/hwglabr2/wiki/Apply-function-to-GRanges-scores-in-genomic-tiles
   forBG<-resGR
   mcols(forBG)<-mcols(forBG)[,c("wormbaseID","score")]
   colnames(mcols(forBG))<-c("name","score")
   seqinfo(forBG)<-ce11seqinfo
   export(forBG,paste0(outPath,"/tracks/",fileNamePrefix,contrastNames[[grp]],
                        "_wt_lfc.bedGraph"),
          format="bedGraph")


   forBW<-disjoin(forBG,ignore.strand=T)
   oldf<-as.data.frame(findOverlaps(forBW,forBG,ignore.strand=T))
   oldf$scorePerSubBp<-forBG$score[oldf$subjectHits]/width(forBG)[oldf$subjectHits]
   oldf$scorePerQuery<-width(forBW)[oldf$queryHits]*oldf$scorePerSubBp
   score<-oldf %>% group_by(queryHits) %>% summarise(score=mean(scorePerQuery))
   forBW$score<-score$score
   export(forBW,paste0(outPath,"/tracks/",fileNamePrefix,contrastNames[[grp]],
                       "_wt_lfc.bw"),
          format="bigwig")

   #######-
   # bed file for significant genes ------------------------------------------
   #######-

   idx<-which(resGR$padj<padjVal)
   forBed<-resGR[idx]
   mcols(forBed)<-mcols(forBed)[,c("wormbaseID","score")]
   colnames(mcols(forBed))<-c("name","score")
   seqinfo(forBed)<-ce11seqinfo
   #NaIdx<-is.na(forBed$score)
   #forBed$score[NaIdx]<-0
   export(forBed,paste0(outPath,"/tracks/",fileNamePrefix,contrastNames[[grp]],
                        "_wt_lfc_p",gsub("^0.","",padjVal),".bedGraph"),
          format="bedGraph")


   export(forBed,paste0(outPath,"/tracks/",fileNamePrefix,contrastNames[[grp]],
                        "_wt_lfc_p",gsub("^0.","",padjVal),".bed"),
          format="bed")





   # MAplots -----------------------------------------------------------------
   ###########-
   # MAplot ALL genes
   ############-

   pdf(file=paste0(outPath,"/plots/",fileNamePrefix,contrastNames[[grp]],
                   "_MAplots_results.pdf"), width=5,height=5,paper="a4")


   plotMA(res, main=paste0(grp," uncorrected LFC, threshold=", padjVal), ylim=c(-3,3), alpha=padjVal)
   plotMA(resLFC, main=paste0(grp," apeglm shrunk LFC, threshold=", padjVal), ylim=c(-3,3), alpha=padjVal)
   #plotCounts(dds, gene=which.min(res$padj), intgroup="sampleType")

   #############-
   # MAplot X chr genes
   #############-

   #chrXgenes<-mcols(dds)$gene[mcols(dds)$chr=="chrX"]
   chrXgenes<-resLFC$wormbaseID[resLFC$chr=="chrX"]
   chrXres<-resLFC[rownames(resLFC) %in% chrXgenes,]
   chrXres05<-chrXres[chrXres$padj<padjVal,]

   if(length(chrXgenes)>0) {
      upOnX<-chrXres05[chrXres05$log2FoldChange>0,]
      write.table(rownames(upOnX), file=paste0(outPath,"/txt/",
                                               fileNamePrefix, contrastNames[[grp]],
                                            "_upOnX_p",padjVal,".csv"),
               row.names=FALSE,col.names=FALSE)

      plotMA(chrXres,main=paste0(grp, " chrX genes, threshold= ", padjVal),
             ylim=c(-4,4),alpha=padjVal)
   }


   #############-
   # MAplotautosomal genes
   #############-
   autosomalGenes<-resLFC$wormbaseID[resLFC$chr %in% c("chrI","chrII","chrIII","chrIV","chrV")]
   #autosomalGenes<-mcols(dds)$gene[mcols(dds)$chr %in% c("chrI","chrII","chrIII","chrIV","chrV")]
   autosomalRes<-resLFC[rownames(resLFC) %in% autosomalGenes,]

   autosomalRes05<- autosomalRes[autosomalRes$padj<padjVal,]

   plotMA(autosomalRes, main=paste0(grp, " autosomal genes, threshold=",
                                    padjVal),ylim=c(-4,4),alpha=padjVal)

   dev.off()


   # Fisher tests ------------------------------------------------------------
   #############-
   # Fisher test of number of up and down genes on X v autosomes
   #############-

   sink(file=paste0(outPath,"/txt/",fileNamePrefix,contrastNames[[grp]],
                    "_logfile.txt"),append=TRUE, type="output")
   upVdownXvA<-matrix(data=c(sum(chrXres05$log2FoldChange>0),
                             sum(chrXres05$log2FoldChange<0),
                             sum(autosomalRes05$log2FoldChange>0),
                             sum(autosomalRes05$log2FoldChange<0)),nrow=2,
                      dimnames=list(group=c("Up","Down"),
                                    chr=c("chrX","chrA")))

   cat("\nFisher Test, up v down:\n")
   print(upVdownXvA)
   print(fisher.test(upVdownXvA))


   #############-
   # Fisher test of number of differentially expressed genes on X v autosomes
   #############-

   testEnrich<-matrix(c(dim(chrXres)[1],dim(chrXres05)[1],
                        dim(autosomalRes)[1],
                        dim(autosomalRes05)[1]),
                      nrow=2,dimnames=list(group=c("NumTotal","NumSig"),chr=c("chrX","chrA")))
   cat("\nFisher Test, enrichment of differentially expressed genes:\n")
   print(testEnrich)
   print(fisher.test(testEnrich))
   sink()


   # boxplots X vs autosomes -------------------------------------------------
   #############-
   # Box plot by X v autosomes
   #############-
   pdf(file=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                   "_boxPlots_expnByChrType.pdf"), width=5,height=5,paper="a4")

   idx<-resLFC$log2FoldChange!=0
   chrType<-factor(rownames(resLFC) %in% chrXgenes)
   levels(chrType)<-c("Autosomal","X chr")
   geneCounts<-table(chrType)

   boxplot(log2FoldChange~chrType, data=resLFC, varwidth=TRUE, outline=FALSE, notch=TRUE,
           main=paste0("Expression changes ", grp), col="grey", ylab="Log2 Fold Change",
           xlab="chromosome type (number of genes)", names=paste(names(geneCounts)," \n(",geneCounts,")",sep=""))
   #stripchart(log2FoldChange~chrType,data=res,method="jitter",vertical=TRUE,pch=20,col="#11115511",cex=0.5,add=TRUE)
   abline(h=0,lty=2,col="blue")

   dev.off()


   #############-
   # Box plot by chromosome
   #############-
   pdf(file=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                   "_boxPlots_expnByChr.pdf"), width=8,height=5,paper="a4")
   chrName<-factor(resLFC$chr)
   geneCounts<-table(chrName)

   boxplot(log2FoldChange~chrName,data=resLFC,varwidth=TRUE,outline=FALSE,notch=TRUE,
           main=paste0("Expression changes ", grp), ylab="log2 Fold Change",
           col=c(rep("grey",5),"purple"),xlab="chromosome (number of genes)",
           names=paste(names(geneCounts)," \n(",geneCounts,")",sep=""))
   #stripchart(log2FoldChange~chrType,data=res,method="jitter",vertical=TRUE,pch=20,col="#11115511",cex=0.5,add=TRUE)
   abline(h=0,lty=2,col="blue")


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
   keyvals[which(resLFC$padj<padjVal & abs(resLFC$log2FoldChange)>lfcVal)]<-'red'
   names(keyvals)[which(resLFC$padj<padjVal & abs(resLFC$log2FoldChange)>lfcVal)]<-paste0('p<',padjVal,' |lfc|>',lfcVal)
   sigUp<-sum(resLFC$padj<padjVal & resLFC$log2FoldChange>lfcVal)
   sigDown<-sum(resLFC$padj<padjVal & resLFC$log2FoldChange< -lfcVal)
   p1<-EnhancedVolcano(resLFC,
                   lab=rownames(resLFC),
                   labSize=0.5,
                   labCol="#11111100",
                   x="log2FoldChange",
                   y="padj",
                   selectLab=rownames(resLFC)[12366],
                   xlim=c(-5.5,5.5),
                   ylim=c(0,65),
                   title= grp,
                   titleLabSize = 16,
                   subtitle=NULL,
                   caption = paste0(sum(!is.na(resLFC$padj)), ' total genes. ',sigUp,
                                    " up, ",sigDown," down."),
                   captionLabSize = 12,
                   pCutoff=padjVal,
                   FCcutoff=lfcVal,
                   xlab=bquote(~Log[2]~'FC'~.(grp)),
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
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                             "_volcanoPlot_allGenes.pdf"), plot=p1,
             device="pdf",path=outPath, width=12,height=12,units="cm")
   } else {
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                          "_volcanoPlot_allGenes.png"), plot=p1,
          device="png",path=outPath, width=12,height=12,units="cm")
   }

   resByChr<-resLFC[order(resLFC$chr),]
   # create custom key-value pairs for 'low', 'chrX', 'autosome' expression by fold-change
   # set the base colour as 'black'
   keyvals <- rep('black', nrow(resByChr))
   # set the base name/label as 'NS'
   names(keyvals) <- rep('NS', nrow(resByChr))
   # modify keyvals for variables with fold change > 2.5
   keyvals[which(resByChr$chr=="chrX")] <- 'red2'
   names(keyvals)[which(resByChr$chr=="chrX")] <- 'chrX'

   # modify keyvals for variables with fold change < -2.5
   keyvals[which(resByChr$chr!="chrX")] <- 'royalblue'
   names(keyvals)[which(resByChr$chr!="chrX")] <- 'autosomes'

   # modify keyvals for variables with fold change < -2.5
   #keyvals[which(abs(resByChr$log2FoldChange)<1 | resByChr$padj>5*10e-3)] <- 'grey10'
   #names(keyvals)[which(abs(resByChr$log2FoldChange)<1 | resByChr$padj>5*10e-3)] <- 'NS'


   #pdf(file=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
   #                "_volcanoPlot_expnByChr.pdf"), width=8,height=6,paper="a4")
   sigUp<-sum(resByChr$padj<padjVal & resByChr$log2FoldChange>lfcVal)
   sigDown<-sum(resByChr$padj<padjVal & resByChr$log2FoldChange< -lfcVal)
   p2<-EnhancedVolcano(resByChr,
                   lab=rownames(resByChr),
                   labSize=0.5,
                   labCol="#11111100",
                   x="log2FoldChange",y="padj",
                   selectLab=rownames(resByChr)[12366],
                   xlim=c(-5.5,5.5),
                   ylim=c(0,65),
                   title= grp,
                   titleLabSize = 16,
                   subtitle=NULL,
                   caption = paste0(sum(!is.na(resLFC$padj)), ' tested genes. ',sigUp, " up, ",sigDown," down."),
                   captionLabSize = 12,
                   pCutoff=padjVal,
                   FCcutoff=lfcVal,
                   xlab=bquote(~Log[2]~'FC'~.(grp)),
                   ylab=bquote(~-Log[10]~adjusted~italic(P)),
                   legendPosition = 'top',
                   legendLabSize = 12,
                   legendIconSize = 3.0,
                   axisLabSize=14,
                   colCustom=keyvals,
                   colAlpha=0.5,
                   pointSize = 1.0)
   #dev.off()
   if(plotPDFs==T){
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                             "_volcanoPlot_autVchrX.pdf"), plot=p2,
             device="pdf",path=outPath,width=12,height=12,units="cm")
   } else {
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                          "_volcanoPlot_autVchrX.png"), plot=p2,
          device="png",path=outPath,width=12,height=12,units="cm")
   }

   if(length(chrXgenes)>0) {
      idx<-resByChr$chr=="chrX"
      #pdf(file=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
      #                "_volcanoPlot_chrX.pdf"), width=8,height=6,paper="a4")
      sigUp<-sum(resByChr$padj[idx]<padjVal & resByChr$log2FoldChange[idx]>lfcVal)
      sigDown<-sum(resByChr$padj[idx]<padjVal & resByChr$log2FoldChange[idx]< -lfcVal)
      p3<-EnhancedVolcano(resByChr[idx,],
                          lab=rownames(resByChr[idx,]),
                          x="log2FoldChange",y="padj",
                          selectLab=rownames(resByChr)[12366],
                          xlim=c(-5.5,5.5),
                          ylim=c(0,65),
                          title= paste0(grp,": chrX genes"),
                          titleLabSize = 16,
                          subtitle=NULL,
                          caption = paste0(sum(idx), ' total genes. ',sigUp, " up, ",sigDown," down."),
                          captionLabSize = 12,
                          pCutoff=padjVal,
                          FCcutoff=lfcVal,
                          xlab=bquote(~Log[2]~'FC'~.(grp)),
                          ylab=bquote(~-Log[10]~adjusted~italic(P)),
                          legendPosition = 'top',
                          legendLabSize = 12,
                          legendIconSize = 3.0,
                          axisLabSize=14,
                          colCustom=keyvals[idx],
                          colAlpha=0.5,
                          pointSize=1.0)
      #dev.off()
      if(plotPDFs==T){
         ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                                "_volcanoPlot_chrX.pdf"), plot=p3,
                device="pdf",path=outPath,width=12,height=12,units="cm")
      } else {
         ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                                "_volcanoPlot_chrX.png"), plot=p3,
                device="png",path=outPath,width=12,height=12,units="cm")
      }
   }

   idx<-resByChr$chr!="chrX"
   #pdf(file=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
   #                "_volcanoPlot_autosomes.pdf"), width=8,height=6,paper="a4")
   sigUp<-sum(resByChr$padj[idx]<padjVal & resByChr$log2FoldChange[idx]>lfcVal)
   sigDown<-sum(resByChr$padj[idx]<padjVal & resByChr$log2FoldChange[idx]< -lfcVal)
   p4<-EnhancedVolcano(resByChr[idx,],
                   lab=rownames(resByChr[idx,]),
                   labSize=0.5,
                   labCol="#11111100",
                   x="log2FoldChange",y="padj",
                   selectLab=rownames(resByChr)[12366],
                   xlim=c(-5.5,5.5),
                   ylim=c(0,65),
                   title= paste0(grp,": autosomal genes"),
                   titleLabSize = 16,
                   subtitle=NULL,
                   caption = paste0(sum(idx), ' total genes. ',sigUp, " up, ",sigDown," down."),
                   captionLabSize = 12,
                   pCutoff=padjVal,
                   FCcutoff=lfcVal,
                   xlab=bquote(~Log[2]~'FC'~.(grp)),
                   ylab=bquote(~-Log[10]~adjusted~italic(P)),
                   legendPosition = 'top',
                   legendLabSize = 12,
                   legendIconSize = 3.0,
                   axisLabSize=14,
                   colCustom=keyvals[idx],
                   colAlpha=0.5,
                   pointSize=1.0)
   #dev.off()
   if (plotPDFs==T) {
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                             "_volcanoPlot_autosomes.pdf"), plot=p4,
             device="pdf",path=outPath,width=12,height=12,units="cm")
   } else {
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                      "_volcanoPlot_autosomes.png"), plot=p4,
          device="png",path=outPath,width=12,height=12,units="cm")
   }



   sink(file=paste0(outPath,"/txt/", fileNamePrefix, contrastNames[[grp]],
                    "_logfile.txt"),append=TRUE, type="output")
   cat("Summary by Chr: \n")
   cat("\np=0.05, LFC=0: \n")
   print(summaryByChr(resLFC,padj=0.05,lfc=0))
   cat("\np=0.05, LFC=0.5: \n")
   print(summaryByChr(resLFC,padj=0.05,lfc=0.5))
   cat("\np=0.05, LFC=1: \n")
   print(summaryByChr(resLFC,padj=0.05,lfc=1))

   cat("\np=0.01, LFC=0: \n")
   print(summaryByChr(resLFC,padj=0.01,lfc=0))
   cat("\np=0.01, LFC=0.5: \n")
   print(summaryByChr(resLFC,padj=0.01,lfc=0.5))
   cat("\np=0.01, LFC=1: \n")
   print(summaryByChr(resLFC,padj=0.01,lfc=1))

   sink()

   if(length(chrXgenes)>0) {
      salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],
                             "_DESeq2_fullResults_p",padjVal,".rds"))
      salmondc<-filterResults(salmon,padj=padjVal,lfc=lfcVal,"gt","chrX", writeTable=F)
      salmondcgr<-metadata[metadata$wormbaseID %in% salmondc$wormbaseID]
      mcols(salmondcgr)<-cbind(mcols(salmondcgr),
                               salmondc[match(salmondcgr$wormbaseID,
                                              salmondc$wormbaseID),c(1:3)])
      salmondcgr
      saveRDS(salmondcgr,file=paste0(outPath,"/rds/",fileNamePrefix, contrastNames[[grp]],
                                     "_chrXup_lfc", lfcVal,"_p",
                                     padjVal, "_gr.rds"))
   }
}



# Volcano - colour by other datasets --------------------------------------

oscillating<-read.delim(paste0(outPath,"/publicData/oscillatingGenes.tsv"),
                        header=T, stringsAsFactors=F)
latorre<-read.delim(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv")) #3235
hsUp<-readRDS(paste0(outPath,"/publicData/hsUp_garrigues2019.rds"))
hsDown<-readRDS(paste0(outPath,"/publicData/hsDown_garrigues2019.rds"))

amplicons<-readRDS(paste0(outPath,"/otherData/ampliconMaxTSSgr.RDS"))


#grp=groupsOI[3]

for(grp in names(contrastNames)){
   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],
                          "_DESeq2_fullResults_p",padjVal,".rds"))

   #### oscillating genes
   bkgrnd='#99999966'
   keyvals<-rep(bkgrnd, nrow(salmon))
   names(keyvals)<-rep('Other',nrow(salmon))
   idx<-(salmon$wormbaseID %in% oscillating$wormbaseID) & (salmon$wormbaseID %in% latorre$wormbaseID)
   keyvals[idx]<-'blue'
   names(keyvals)[idx]<-"Both"
   idx<-(salmon$wormbaseID %in% oscillating$wormbaseID) & !(salmon$wormbaseID %in% latorre$wormbaseID)
   keyvals[idx]<-'red'
   names(keyvals)[idx]<-"Meeuse(2020)"
   idx<-(salmon$wormbaseID %in% latorre$wormbaseID) & !(salmon$wormbaseID %in% oscillating$wormbaseID)
   keyvals[idx]<-'green'
   names(keyvals)[idx]<-"Latorre(2015)"
   sigUp<-sum(rowSums(cbind(salmon$padj< padjVal,
                            salmon$log2FoldChange>lfcVal,
                            keyvals!=bkgrnd), na.rm=T)==3, na.rm=T)
   sigDown<-sum(rowSums(cbind(salmon$padj< padjVal,
                              salmon$log2FoldChange< -lfcVal,
                              keyvals!=bkgrnd), na.rm=T)==3, na.rm=T)
   p1<-EnhancedVolcano(salmon,
                       lab=salmon$publicID,
                       labSize=0.5,
                       labCol="#11111100",
                       x="log2FoldChange",
                       y="padj",
                       selectLab=salmon$publicID[12366],
                       xlim=c(-5.5,5.5),
                       ylim=c(0,65),
                       title= grp,
                       subtitle=NULL,
                       caption = paste0(sum(keyvals!=bkgrnd), ' oscillating genes. ',sigUp, " up, ",sigDown," down."),
                       captionLabSize = 12,
                       pCutoff=padjVal,
                       FCcutoff=lfcVal,
                       xlab=bquote(~Log[2]~'fold change'~.(grp)),
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
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                             "_volcanoPlot_oscillatingGenes.pdf"), plot=p1,
             device="pdf",path=outPath, width=12,height=12,units="cm")
   } else {
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                             "_volcanoPlot_oscillatingGenes.png"), plot=p1,
             device="png",path=outPath, width=12,height=12,units="cm")
   }


   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],
                          "_DESeq2_fullResults_p",padjVal,".rds"))



   #### heat shock genes
   myCols<-c("#11111100","red","red") # background, dataset1, dataset2
   keyvals<-rep(myCols[1], nrow(salmon))
   names(keyvals)<-rep('Other',nrow(salmon))
   idx<-salmon$wormbaseID %in% hsUp$wormbaseID | salmon$wormbaseID %in% hsDown$wormbaseID
   keyvals[idx]<-myCols[2]
   names(keyvals)[idx]<-"Heatshock"
   sigUp<-sum(rowSums(cbind(salmon$padj< padjVal,
                            salmon$log2FoldChange>lfcVal,
                            keyvals==myCols[2]), na.rm=T)==3, na.rm=T)
   sigDown<-sum(rowSums(cbind(salmon$padj< padjVal,
                              salmon$log2FoldChange< -lfcVal,
                              keyvals==myCols[2]), na.rm=T)==3, na.rm=T)
   p1<-EnhancedVolcano(salmon,
                       lab=salmon$publicID,
                       labSize=0.5,
                       labCol=myCols[1],
                       x="log2FoldChange",
                       y="padj",
                       selectLab=salmon$publicID[12366],
                       xlim=c(-5.5,5.5),
                       ylim=c(0,65),
                       title= grp,
                       subtitle=NULL,
                       caption = paste0(sum(keyvals==myCols[2]), ' heatshock genes (Garrigues 2019). ',sigUp, " up, ",sigDown," down."),
                       captionLabSize = 12,
                       pCutoff=padjVal,
                       FCcutoff=lfcVal,
                       xlab=bquote(~Log[2]~'fold change'~.(grp)),
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
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                             "_volcanoPlot_hsGenes.pdf"), plot=p1,
             device="pdf",path=outPath, width=12,height=12,units="cm")
   } else {
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                             "_volcanoPlot_hsGenes.png"), plot=p1,
             device="png",path=outPath, width=12,height=12,units="cm")
   }


   #### amplicon genes
   amp<-salmon[salmon$wormbaseID %in% amplicons$WBgeneID,]
   keyvals<-rep('black', nrow(amp))
   names(keyvals)<-rep('Other',nrow(amp))
   idxX<-amp$chr=="chrX"
   idxA<-amp$chr!="chrX"
   keyvals[idxX]<-'red'
   keyvals[idxA]<-'blue'
   names(keyvals)[idxX]<-"chrX"
   names(keyvals)[idxA]<-"Autosomes"
   sigUp<-sum(rowSums(cbind(amp$padj< padjVal,
                            amp$log2FoldChange>lfcVal), na.rm=T)==2, na.rm=T)
   sigDown<-sum(rowSums(cbind(amp$padj< padjVal,
                              amp$log2FoldChange< -lfcVal), na.rm=T)==2, na.rm=T)
   p1<-EnhancedVolcano(amp,
                       lab=amp$publicID,
                       labSize=2.5,
                       labCol="#11111144",
                       x="log2FoldChange",
                       y="padj",
                       selectLab=amp$publicID[amp$padj< padjVal &
                                                 abs(amp$log2FoldChange)>lfcVal],
                       drawConnectors = TRUE,
                       widthConnectors = 0.2,
                       colConnectors = "#11111133",
                       lengthConnectors = unit(0.01,'snpc'),
                       xlim=c(-5.5,5.5),
                       ylim=c(0,65),
                       title= grp,
                       subtitle=NULL,
                       caption = paste0(sum(keyvals!="black"), ' amplicon genes: ',sigUp, " up, ",sigDown," down."),
                       captionLabSize = 12,
                       pCutoff=padjVal,
                       FCcutoff=lfcVal,
                       xlab=bquote(~Log[2]~'fold change'~.(grp)),
                       ylab=bquote(~-Log[10]~adjusted~italic(P)),
                       legendPosition = 'top',
                       legendLabSize = 12,
                       legendIconSize = 3.0,
                       axisLabSize=14,
                       colCustom=keyvals,
                       colAlpha=0.5,
                       pointSize = 1.0)

   if(plotPDFs==T){
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                             "_volcanoPlot_ampliconGenes.pdf"), plot=p1,
             device="pdf",path=outPath, width=12,height=12,units="cm")
   } else {
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                             "_volcanoPlot_ampliconGenes.png"), plot=p1,
             device="png",path=outPath, width=12,height=12,units="cm")
   }

}


##########-
# LFC density plots------
##########-

#padjVals=c(0.05,0.99)
lfcDensity<-NULL
for(grp in names(contrastNames)) {
   #######-
   # vary filtering threshold ----------------------------------------
   #######-
   resLFC<-readRDS(paste0(outPath,"/rds/", fileNamePrefix, contrastNames[[grp]],
                          "_DESeq2_fullResults_p",padjVal,".rds"))
   ## Count significant genes at different thresholds and plot
   thresholds<-varyThreshold1(resLFC,
                             pval=padjVal, lfcVals=c(0,0.25,0.5,1),
                             direction=c("both","lt","gt"),
                             chr=c("all","chrX","autosomes"))

   thresholds<-thresholds[order(thresholds$group,thresholds$padj,thresholds$direction,thresholds$lfc),]
   write.csv(thresholds,file=paste0(outPath,"/txt/",fileNamePrefix, contrastNames[[grp]],
                                    "_numSignificant_p",padjVal,".csv"),quote=F,row.names=F)

   ####### plot percentSignificant
   p1<-ggplot(data=thresholds,aes(x=as.factor(lfc),y=percentSignificant))+
      facet_grid(rows=vars(direction),cols=vars(chr))+
      ggtitle(grp)+geom_bar(stat="identity")+
      xlab("log2 fold change threshold")+ylab("Percent significant genes") +
      geom_text(aes(label=paste0(round(percentSignificant,0),"%\n",numSignificant)),
                vjust=0,size=3,lineheight=0.9)+
      ylim(0,1.2*max(thresholds$percentSignificant))

   ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                          "_thresholds_percentSig_p",padjVal,".png"), plot=p1,
          device="png",path=outPath,width=14,height=12,units="cm")

   p2<-ggplot(data=thresholds,aes(x=as.factor(lfc),y=percentSigGt10))+
      facet_grid(rows=vars(direction),cols=vars(chr))+
      ggtitle(grp)+geom_bar(stat="identity")+
      xlab("log2 fold change threshold")+ylab("Percent significant genes") +
      geom_text(aes(label=paste0(round(percentSigGt10,0),"%\n",numSigGt10)),
                vjust=0,size=3,lineheight=0.9)+
      ylim(0,1.2*max(thresholds$percentSigGt10))

   ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, contrastNames[[grp]],
                          "_thresholds_percentSigGt10_p",padjVal,".png"), plot=p2,
          device="png",path=outPath,width=14,height=12,units="cm")


   ### plot distribution in full-----
   dd<-getDensity1(resLFC, pval=padjVal, breaks=c(seq(0,2,0.05),Inf),
                  chr="all", asCounts=F)[[1]]
   dd$group<-grp
   if(is.null(lfcDensity)){
      lfcDensity<-dd
   } else {
      lfcDensity<-rbind(lfcDensity,dd)
   }


   if(length(chrXgenes)>0){
      dd<-getDensity1(resLFC, pval=padjVal,
                     breaks=c(seq(0,2,0.05),Inf), chr="chrX",
                     direction="gt", asCounts=F)
      rawX<-dd[[2]]
      dd<-dd[[1]]
      dd$group<-paste0(grp,"_chrX")
      lfcDensity<-rbind(lfcDensity,dd)
   }

   dd<-getDensity1(resLFC, pval=padjVal,
                  breaks=c(seq(0,2,0.05),Inf), chr="autosomes",
                  direction="both", asCounts=F)[[1]]
   dd$group<-paste0(grp,"_chrA")
   lfcDensity<-rbind(lfcDensity,dd)
}

lfcDensity$breaks<-gsub(",","-",gsub("\\(|\\]","",lfcDensity$breaks))
lfcDensity$breaks<-factor(lfcDensity$breaks)

densityGroups<-length(unique(lfcDensity$group))
if(densityGroups>9){
   lfcDensity<-lfcDensity[grep("_chr[A|X]",lfcDensity$group),]
   densityGroups<-length(unique(lfcDensity$group))
}

p<-ggplot(data=lfcDensity,aes(x=breaks,y=counts)) + facet_grid(rows=vars(group),cols=vars(pvals),switch = "y")+
   geom_bar(stat="identity") +theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                     strip.text.y.left = element_text(angle = 0,hjust=0)) +
   xlab("Absolute log2 fold change bins")+ylab("Density")


p1<-p+geom_vline(aes(xintercept = 10.5),color="red")+
   annotate("text",label="lfc=0.5",size=3,x=11,y=0.95*max(lfcDensity$counts),
            hjust=0.1,color="red")

p2<-p+geom_vline(aes(xintercept = 5.5),color="red")+
   annotate("text",label="lfc=0.25",size=3,x=5.5,y=0.95*max(lfcDensity$counts),
            hjust=-0.1,color="red")

if(plotPDFs==T){
   ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix,
                          "lfcValueDistribution_p",padjVal,"_0.5.pdf"), plot=p1,
          device="pdf",path=outPath, width=19,height=29*densityGroups/18,units="cm")
   ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix,
                          "lfcValueDistribution_p",padjVal,"_0.25.pdf"), plot=p2,
          device="pdf",path=outPath, width=19,height=29*densityGroups/18,units="cm")
} else {
   ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix,
                          "lfcValueDistribution_p",padjVal,"_0.5.png"), plot=p1,
          device="png",path=outPath, width=19,height=29*densityGroups/18,units="cm")
   ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix,
                          "lfcValueDistribution_p",padjVal,"_0.25.png"), plot=p2,
          device="png",path=outPath, width=19,height=29*densityGroups/18,units="cm")
}




########################-
## ECDF of data -----
########################-

# using automatic filtering threshold
sigTables<-list()
localPadj=padjVal
localLFC=0
for (grp in names(contrastNames)){
   print(grp)
   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],
                          "_DESeq2_fullResults_p",padjVal,".rds"))
   salmon<-salmon[!is.na(salmon$padj),]
   #nrow(filterResults(salmon,padj=0.05,lfc=0.5,direction="lt",chr="autosomes"))
   print(paste0(nrow(salmon)," genes before filtering"))
   print(paste0(sum(is.na(salmon$log2FoldChange))," have log2FoldChange that is NA"))
   #salmon$expressed<-sum(salmon$baseMean>10)
   sigTables[[grp]]<-as.data.frame(salmon) #[salmon$baseMean>10,]
   print(paste0(nrow(sigTables[[grp]])," genes after automatic threshold filter"))
}

SMC<-rep(names(sigTables),lapply(sigTables,nrow))
sig<-do.call(rbind,sigTables)
sig$SMC<-SMC
#sum(is.na(sig$padj))
#sig<-sig[!is.na(sig$padj),]
#sig$SMC<-factor(SMC)
table(sig$SMC)
sig$XvA<-"Autosomes"
sig$XvA[sig$chr=="chrX"]<-"chrX"
#sig$XvA<-factor(sig$XvA)
table(sig$XvA)
sig$upVdown<-"0"
sig$upVdown[sig$log2FoldChange<0]<-"down"
sig$upVdown[sig$log2FoldChange>0]<-"up"
#sig$upVdown<-factor(sig$upVdown,levels=c("0","up","down"))
table(sig$upVdown)
row.names(sig)<-NULL
SMC<-NULL

# check if datasets have chrX genes included
includeChrX<-"chrX" %in% unlist(lapply(sigTables,"[","chr"))
options(tibble.width=Inf)
dd1<-sig %>% dplyr::filter(padj<localPadj) %>%
   dplyr::group_by(SMC,upVdown,XvA) %>%
   dplyr::mutate(ecd=ecdf(abs(log2FoldChange))(abs(log2FoldChange)))

ss1<-sig %>% dplyr::group_by(SMC,XvA) %>% dplyr::mutate(expOnChr=n()) %>%
   dplyr::filter(padj<localPadj) %>%
   dplyr::mutate(sigOnChr=n())%>% dplyr::group_by(SMC,XvA,upVdown) %>%
   dplyr::mutate(sigInGrp=n(), ecd=ecdf(abs(log2FoldChange))(abs(log2FoldChange))) %>%
   dplyr::summarize(expOnChr=unique(expOnChr),
             sigOnChr=unique(sigOnChr),
             sigInGrp=unique(sigInGrp),
             qnt0=1-ecdf(abs(log2FoldChange))(c(0)),
             qnt25=1-ecdf(abs(log2FoldChange))(0.25),
             qnt50=1-ecdf(abs(log2FoldChange))(0.5),
             fractionSig=sigInGrp/sigOnChr,
             count0Sig=sigInGrp*qnt0,
             countq25Sig=sigInGrp*qnt25,
             countq50Sig=sigInGrp*qnt50,
             percentq0Exp=sigInGrp*qnt0*100/expOnChr,
             percentq25Exp=sigInGrp*qnt25*100/expOnChr,
             percentq50Exp=sigInGrp*qnt50*100/expOnChr,
             .groups="keep")
write.table(ss1,file=paste0(outPath,"/txt/ecdf_lfcThresholds_p",padjVal,".tsv"),
            sep="\t",row.names=F,quote=F)


p<-ggplot(dd1, aes(x=abs(log2FoldChange),y=ecd,color=SMC,linetype=XvA)) +
   geom_line(size=0.9)+ facet_wrap(vars(upVdown),nrow=2)+
   theme_classic() + xlim(c(0,1.5)) +
   xlab("Absolute log2 fold change")+ylab("Fraction significant genes rejected")
#p

#stat_ecdf(aes(colour=SMC,linetype=XvA),alpha=0.7)
p1<-p+geom_vline(aes(xintercept = 0.5), color="grey") +
   annotate("text",label="0.5",size=3, x=0.5, y=0,hjust=-0.05,color="grey") +
   geom_vline(aes(xintercept = 0.25), color="grey") +
   annotate("text",label="0.25",size=3, x=0.25, y=0,hjust=-0.05,color="grey")
#p1

if(plotPDFs==T){
   ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix,
                          "lfcValueCDF_p",padjVal,".pdf"), plot=p1,
          device="pdf",path=outPath, width=10,height=10,units="cm")
} else {
   ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix,
                          "lfcValueCDF_p",padjVal,".png"), plot=p1,
          device="png",path=outPath, width=10,height=10,units="cm")
}





# using baseMean>10 threshold
sigTables<-list()
localPadj=padjVal
localLFC=0
for (grp in names(contrastNames)){
   print(grp)
   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],
                          "_DESeq2_fullResults_p",padjVal,".rds"))
   salmon<-salmon[!is.na(salmon$padj),]
   #nrow(filterResults(salmon,padj=0.05,lfc=0.5,direction="lt",chr="autosomes"))
   print(paste0(nrow(salmon)," genes before filtering"))
   print(paste0(sum(is.na(salmon$log2FoldChange))," have log2FoldChange that is NA"))
   #salmon$expressed<-sum(salmon$baseMean>10)
   sigTables[[grp]]<-as.data.frame(salmon[salmon$baseMean>10,])
   print(paste0(nrow(sigTables[[grp]])," genes after >10 baseMean filter"))
}

SMC<-rep(names(sigTables),lapply(sigTables,nrow))
sig<-do.call(rbind,sigTables)
sig$SMC<-SMC
#sum(is.na(sig$padj))
#sig<-sig[!is.na(sig$padj),]
#sig$SMC<-factor(SMC)
table(sig$SMC)
sig$XvA<-"Autosomes"
sig$XvA[sig$chr=="chrX"]<-"chrX"
#sig$XvA<-factor(sig$XvA)
table(sig$XvA)
sig$upVdown<-"0"
sig$upVdown[sig$log2FoldChange<0]<-"down"
sig$upVdown[sig$log2FoldChange>0]<-"up"
#sig$upVdown<-factor(sig$upVdown,levels=c("0","up","down"))
table(sig$upVdown)
row.names(sig)<-NULL
SMC<-NULL

# check if datasets have chrX genes included
includeChrX<-"chrX" %in% unlist(lapply(sigTables,"[","chr"))
options(tibble.width=Inf)
dd1<-sig %>% filter(padj<localPadj) %>%
   dplyr::group_by(SMC,upVdown,XvA) %>%
   dplyr::mutate(ecd=ecdf(abs(log2FoldChange))(abs(log2FoldChange)))

ss1<-sig %>% dplyr::group_by(SMC,XvA) %>% dplyr::mutate(expOnChr=n()) %>%
   dplyr::filter(padj<localPadj) %>%
   dplyr::mutate(sigOnChr=n())%>% dplyr::group_by(SMC,XvA,upVdown) %>%
   dplyr::mutate(sigInGrp=n(), ecd=ecdf(abs(log2FoldChange))(abs(log2FoldChange))) %>%
   dplyr::summarize(expOnChr=unique(expOnChr),
             sigOnChr=unique(sigOnChr),
             sigInGrp=unique(sigInGrp),
             qnt0=1-ecdf(abs(log2FoldChange))(c(0)),
             qnt25=1-ecdf(abs(log2FoldChange))(0.25),
             qnt50=1-ecdf(abs(log2FoldChange))(0.5),
             fractionSig=sigInGrp/sigOnChr,
             count0Sig=sigInGrp*qnt0,
             countq25Sig=sigInGrp*qnt25,
             countq50Sig=sigInGrp*qnt50,
             percentq0Exp=sigInGrp*qnt0*100/expOnChr,
             percentq25Exp=sigInGrp*qnt25*100/expOnChr,
             percentq50Exp=sigInGrp*qnt50*100/expOnChr,
             .groups="keep")
write.table(ss1,file=paste0(outPath,"/txt/ecdf_lfcThresholds_p",padjVal,"_gt10.tsv"),
            sep="\t",row.names=F,quote=F)


p<-ggplot(dd1, aes(x=abs(log2FoldChange),y=ecd,color=SMC,linetype=XvA)) +
   geom_line(size=0.9)+ facet_wrap(vars(upVdown),nrow=2)+
   theme_classic() + xlim(c(0,1.5)) +
   xlab("Absolute log2 fold change")+ylab("Fraction significant genes rejected")
#p

#stat_ecdf(aes(colour=SMC,linetype=XvA),alpha=0.7)
p1<-p+geom_vline(aes(xintercept = 0.5), color="grey") +
   annotate("text",label="0.5",size=3, x=0.5, y=0,hjust=-0.05,color="grey") +
   geom_vline(aes(xintercept = 0.25), color="grey") +
   annotate("text",label="0.25",size=3, x=0.25, y=0,hjust=-0.05,color="grey")
#p1

if(plotPDFs==T){
   ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix,
                          "lfcValueCDF_p",padjVal,"gt10.pdf"), plot=p1,
          device="pdf",path=outPath, width=10,height=10,units="cm")
} else {
   ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix,
                          "lfcValueCDF_p",padjVal,"gt10.png"), plot=p1,
          device="png",path=outPath, width=10,height=10,units="cm")
}





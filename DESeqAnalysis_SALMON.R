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

#library(REBayes)
# get funciton for converting gene names from WBID to publicID
source("~/Documents/MeisterLab/GenomeVer/geneNameConversion/convertingGeneNamesFunction1.R")
source("functions.R")
####
### some variables
####
plotPDFs=F
fileNamePrefix="salmon_"
filterPrefix="noOsc_"
filterData=T
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

fileNames<-paste0(outPath,"/salmon/mRNA/",sampleNames,"/quant.sf")

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


# Create metadata object --------------------------------------------------
###############################################################-
### create metadata
###############################################################-

if(!file.exists(paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.",
                       genomeVer, ".annotations.sqlite"))){
  dir.create(paste0(genomeDir,"/annotations"),recursive=T)
  system(paste0("curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/gff/c_elegans.PRJNA13758.",
                genomeVer,".annotations.gff3.gz -o ",genomeDir,
                "/annotations/c_elegans.PRJNA13758.",genomeVer,
                ".annotations.gff3.gz"))

  system(paste0("gunzip ",genomeDir,"/annotations/c_elegans.PRJNA13758.",
                genomeVer,".annotations.gff3.gz"))
  si<-seqinfo(Celegans)
  genome(si)<-genomeVer
  seqnames(si)<-gsub("M","MtDNA",gsub("chr","",seqnames(si)))
  wstxdb<-makeTxDbFromGFF(file=paste0(genomeDir,
                                      "/annotations/c_elegans.PRJNA13758.",
                                      genomeVer,".annotations.gff3"),
                          format="gff3",organism="Caenorhabditis elegans",
                          chrominfo=si)

  saveDb(wstxdb,paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.", genomeVer,
                     ".annotations.sqlite"))
  file.remove(paste0(genomeDir,"/annotations/c_elegans.PRJNA13758.",
                     genomeVer, ".annotations.gff3"))
}

# load a txdb of wormbase data and create a tx2gene object
txdb<-loadDb(paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.", genomeVer,
                    ".annotations.sqlite"))
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

columns(txdb)
keytypes(txdb)
TxptByGene<-transcriptsBy(txdb, by = "gene")
length(TxptByGene)

geneGR<-unlist(range(TxptByGene))
mcols(geneGR)$wormbase<-names(geneGR)
genedf<-as.data.frame(geneGR)


# download gene id data from simplemine: https://wormbase.org/tools/mine/simplemine.cgi
# for entrez ids, load wormbaseID column into https://david.ncifcrf.gov/conversion.jsp
geneIDs<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/simplemine_WS278_geneID.txt")
david<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/david_wbid2entrez_WS278.txt")

metadata<-inner_join(geneIDs, genedf,by=c("WormBase.Gene.ID"="wormbase")) %>%
  dplyr::select(WormBase.Gene.ID,Public.Name,Sequence.Name,seqnames,start, end, strand) %>%
  collect %>% GenomicRanges::GRanges()

names(mcols(metadata))<-c("wormbaseID","publicID","sequenceID")

i<-which(metadata$wormbaseID %in% david$From)
j<-match(metadata$wormbaseID[i],david$From)
metadata$entrezID<-NA
metadata$entrezID[i]<-david$To[j]

#seqinfo(metadata)<-wbseqinfo
seqlevelsStyle(metadata)<-"ucsc"
seqinfo(metadata)<-ce11seqinfo
metadata<-sort(metadata)

saveRDS(metadata,paste0(outPath,"/wbGeneGR_WS275.rds"))

###############################################################-
# Import into DESeq2 ------------------------------------------------------
###############################################################-



# import the count matrices
txi<-tximport(sampleTable$fileName,type="salmon",tx2gene=tx2gene)

# read samples into DESeq2
dds <- DESeqDataSetFromTximport(txi=txi,
                                colData=sampleTable,
                          design=~replicate+lane+SMC)



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
vsd <- vst(dds, blind=TRUE)
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
                decreasing=TRUE)[1:500]
df <- as.data.frame(colData(dds)[,c("strain","replicate","lane")])
rownames(df)<-colData(dds)$sampleName
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df, main="Top 500 expressed genes - no clustering")

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=TRUE, annotation_col=df, main="Top 500 expressed genes - clustered by columns")



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


pThresh=0.05
LFCthresh=0
#res=list()
#resLFC=list()

####### filter genes
if(filterData){
   oscillating<-read.delim(paste0(outPath,"/oscillatingGenes.tsv"),header=T,
                        stringsAsFactors=F)
   toFilter<-oscillating$WB_ID

   fileNamePrefix=filterPrefix
}

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
   res$wormbaseID<-rownames(res)
   idx<-match(rownames(res),metadata$wormbaseID)
   res$chr<-factor(seqnames(metadata),levels=paste0("chr",c("I","II","III","IV","V","X")))[idx]
   res$start<-as.vector(start(metadata))[idx]
   res$end<-as.vector(end(metadata))[idx]
   res$strand<-as.vector(strand(metadata))[idx]

   # shrink LFC estimates
   #resultsNames(dds) # to get names of coefficients
   resLFC<-lfcShrink(dds,coef=paste0("SMC_",grp,"_vs_",controlGrp), type="apeglm", res=res)
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

   # remove filtered genes
   if(filterData){
      idx<-resLFC$wormbaseID %in% toFilter
      resLFC<-resLFC[!idx,]
   }

   saveRDS(resLFC,file=paste0(outPath,"/rds/", fileNamePrefix, grp,
                              "_DESeq2_fullResults.rds"))

   #export csv with ordered results
   write.csv(resLFC[order(resLFC$padj),],
             file=paste0(outPath,"/txt/", fileNamePrefix,grp,
                         "_DESeq2_resultsTable.csv"),
             quote=F,row.names=F)

   # remove NAs
   res<-na.omit(res)
   resLFC<-na.omit(resLFC)



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
   gene.kept <- rownames(resLFC)[resLFC$padj <= pThresh & !is.na(resLFC$padj) & abs(resLFC$log2FoldChange)>0.5]

   # Retrieve the normalized counts for gene of interest
   countTable.kept <- log2(counts(dds) + epsilon)[gene.kept, ]
   dim(countTable.kept)
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
             main=paste0(grp," changed genes (p<",pThresh,", lfc>0.5)"))

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
              main=g, las=2, cex.names=1,names.arg=colData(dds)$sampleName)
      legend("topleft",legend=levels(colData(dds)$SMC),fill=c(1,2), cex=0.7)
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
   export(forBG,paste0(outPath,"/tracks/",fileNamePrefix,grp,
                        "_wt_lfc.bedGraph"),
          format="bedGraph")


   forBW<-disjoin(forBG,ignore.strand=T)
   oldf<-as.data.frame(findOverlaps(forBW,forBG,ignore.strand=T))
   oldf$scorePerSubBp<-forBG$score[oldf$subjectHits]/width(forBG)[oldf$subjectHits]
   oldf$scorePerQuery<-width(forBW)[oldf$queryHits]*oldf$scorePerSubBp
   score<-oldf %>% group_by(queryHits) %>% summarise(score=mean(scorePerQuery))
   forBW$score<-score$score
   export(forBW,paste0(outPath,"/tracks/",fileNamePrefix,grp,
                       "_wt_lfc.bw"),
          format="bigwig")

   #######-
   # bed file for significant genes ------------------------------------------
   #######-

   idx<-which(resGR$padj<0.05)
   forBed<-resGR[idx]
   mcols(forBed)<-mcols(forBed)[,c("wormbaseID","score")]
   colnames(mcols(forBed))<-c("name","score")
   seqinfo(forBed)<-ce11seqinfo
   #NaIdx<-is.na(forBed$score)
   #forBed$score[NaIdx]<-0
   export(forBed,paste0(outPath,"/tracks/",fileNamePrefix,grp,
                        "_wt_lfc_p",gsub("^0.","",pThresh),".bedGraph"),
          format="bedGraph")


   export(forBed,paste0(outPath,"/tracks/",fileNamePrefix,grp,
                        "_wt_lfc_p",gsub("^0.","",pThresh),".bed"),
          format="bed")





   # MAplots -----------------------------------------------------------------
   ###########-
   # MAplot ALL genes
   ############-

   pdf(file=paste0(outPath,"/plots/",fileNamePrefix,grp,
                   "_MAplots_results.pdf"), width=5,height=5,paper="a4")


   plotMA(res, main=paste0(grp," uncorrected LFC, threshold=", pThresh), ylim=c(-3,3), alpha=pThresh)
   plotMA(resLFC, main=paste0(grp," apeglm shrunk LFC, threshold=", pThresh), ylim=c(-3,3), alpha=pThresh)
   #plotCounts(dds, gene=which.min(res$padj), intgroup="sampleType")

   #############-
   # MAplot X chr genes
   #############-

   #chrXgenes<-mcols(dds)$gene[mcols(dds)$chr=="chrX"]
   chrXgenes<-resLFC$wormbaseID[resLFC$chr=="chrX"]
   chrXres<-resLFC[rownames(resLFC) %in% chrXgenes,]

   chrXres05<-chrXres[chrXres$padj<pThresh,]


   upOnX<-chrXres05[chrXres05$log2FoldChange>0,]
   write.table(rownames(upOnX), file=paste0(outPath,"/txt/",fileNamePrefix, grp,
                                            "_upOnX_p",pThresh,".csv"),
               row.names=FALSE,col.names=FALSE)

   plotMA(chrXres,main=paste0(grp, " chrX genes, threshold= ", pThresh),ylim=c(-4,4),alpha=pThresh)



   #############-
   # MAplotautosomal genes
   #############-
   autosomalGenes<-resLFC$wormbaseID[resLFC$chr %in% c("chrI","chrII","chrIII","chrIV","chrV")]
   #autosomalGenes<-mcols(dds)$gene[mcols(dds)$chr %in% c("chrI","chrII","chrIII","chrIV","chrV")]
   autosomalRes<-resLFC[rownames(resLFC) %in% autosomalGenes,]

   autosomalRes05<- autosomalRes[autosomalRes$padj<pThresh,]

   plotMA(autosomalRes, main=paste0(grp, " autosomal genes, threshold=",pThresh),ylim=c(-4,4),alpha=pThresh)

   dev.off()


   # Fisher tests ------------------------------------------------------------
   #############-
   # Fisher test of number of up and down genes on X v autosomes
   #############-

   sink(file=paste0(outPath,"/txt/",fileNamePrefix,grp,
                    "_logfile.txt"),append=TRUE, type="output")
   upVdownXvA<-matrix(data=c(sum(chrXres05$log2FoldChange>0),
                             sum(chrXres05$log2FoldChange<0),
                             sum(autosomalRes05$log2FoldChange>0),
                             sum(autosomalRes05$log2FoldChange<0)),nrow=2,dimnames=list(group=c("Up","Down"),chr=c("chrX","chrA")))

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
   pdf(file=paste0(outPath,"/plots/",fileNamePrefix, grp,
                   "_boxPlots_expnByChrType.pdf"), width=5,height=5,paper="a4")

   chrType<-factor(rownames(resLFC) %in% chrXgenes)
   levels(chrType)<-c("Autosomal","X chr")
   geneCounts<-table(chrType)

   boxplot(log2FoldChange~chrType, data=resLFC, varwidth=TRUE, outline=FALSE, notch=TRUE,
           main=paste0("Expression changes after cleavage of ", grp), col="grey", ylab="Log2 Fold Change",
           xlab="chromosome type (number of genes)", names=paste(names(geneCounts)," \n(",geneCounts,")",sep=""))
   #stripchart(log2FoldChange~chrType,data=res,method="jitter",vertical=TRUE,pch=20,col="#11115511",cex=0.5,add=TRUE)
   abline(h=0,lty=2,col="blue")
   dev.off()


   #############-
   # Box plot by chromosome
   #############-
   pdf(file=paste0(outPath,"/plots/",fileNamePrefix, grp,
                   "_boxPlots_expnByChr.pdf"), width=8,height=5,paper="a4")
   chrName<-factor(resLFC$chr)
   geneCounts<-table(chrName)

   boxplot(log2FoldChange~chrName,data=resLFC,varwidth=TRUE,outline=FALSE,notch=TRUE,
           main=paste0("Expression changes after cleavage of ", grp), ylab="log2 Fold Change",
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
   keyvals[which(resLFC$padj<pThresh & abs(resLFC$log2FoldChange)>1)]<-'red'
   names(keyvals)[which(resLFC$padj<pThresh & abs(resLFC$log2FoldChange)>1)]<-paste0('p<',pThresh,' |lfc|>1')
   sigUp<-sum(resLFC$padj<pThresh & resLFC$log2FoldChange>1)
   sigDown<-sum(resLFC$padj<pThresh & resLFC$log2FoldChange< -1)
   p1<-EnhancedVolcano(resLFC,
                   lab=rownames(resLFC),
                   labSize=0.5,
                   labCol="#11111100",
                   x="log2FoldChange",
                   y="padj",
                   selectLab=rownames(resLFC)[12366],
                   xlim=c(-5.5,5.5),
                   ylim=c(0,65),
                   title= paste0(grp," vs ", controlGrp),
                   subtitle=NULL,
                   caption = paste0(nrow(resLFC), ' expressed genes. ',sigUp, " up, ",sigDown," down."),
                   captionLabSize = 12,
                   pCutoff=pThresh,
                   FCcutoff=1.0,
                   xlab=bquote(~Log[2]~'fold change'~.(grp)~'/'~.(controlGrp)),
                   ylab=bquote(~-Log[10]~adjusted~italic(P)),
                   #.legend=c('NS','P & Log2 FC'),
                   #legendLabels=c('NS', expression(p-value<pThresh~and~log[2]~FC>1)),
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


   #pdf(file=paste0(outPath,"/plots/",fileNamePrefix, grp,
   #                "_volcanoPlot_expnByChr.pdf"), width=8,height=6,paper="a4")
   sigUp<-sum(resByChr$padj<pThresh & resByChr$log2FoldChange>1)
   sigDown<-sum(resByChr$padj<pThresh & resByChr$log2FoldChange< -1)
   p2<-EnhancedVolcano(resByChr,
                   lab=rownames(resByChr),
                   labSize=0.5,
                   labCol="#11111100",
                   x="log2FoldChange",y="padj",
                   selectLab=rownames(resByChr)[12366],
                   xlim=c(-5.5,5.5),
                   ylim=c(0,65),
                   title= paste0(grp," vs ",controlGrp),
                   subtitle=NULL,
                   caption = paste0(nrow(resLFC), ' expressed genes. ',sigUp, " up, ",sigDown," down."),
                   captionLabSize = 12,
                   pCutoff=pThresh,
                   FCcutoff=1.0,
                   xlab=bquote(~Log[2]~'fold change'~.(grp)~'/'~.(controlGrp)),
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
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                             "_volcanoPlot_autVchrX.pdf"), plot=p2,
             device="pdf",path=outPath,width=12,height=12,units="cm")
   } else {
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                          "_volcanoPlot_autVchrX.png"), plot=p2,
          device="png",path=outPath,width=12,height=12,units="cm")
   }
   idx<-resByChr$chr=="chrX"
   #pdf(file=paste0(outPath,"/plots/",fileNamePrefix, grp,
   #                "_volcanoPlot_chrX.pdf"), width=8,height=6,paper="a4")
   sigUp<-sum(resByChr$padj[idx]<pThresh & resByChr$log2FoldChange[idx]>1)
   sigDown<-sum(resByChr$padj[idx]<pThresh & resByChr$log2FoldChange[idx]< -1)
   p3<-EnhancedVolcano(resByChr[idx,],
                   lab=rownames(resByChr[idx,]),
                   x="log2FoldChange",y="padj",
                   selectLab=rownames(resByChr)[12366],
                   xlim=c(-5.5,5.5),
                   ylim=c(0,65),
                   title= paste0(grp," vs ",controlGrp,": chrX genes"),
                   subtitle=NULL,
                   caption = paste0(nrow(resLFC), ' expressed genes. ',sigUp, " up, ",sigDown," down."),
                   captionLabSize = 12,
                   pCutoff=pThresh,
                   FCcutoff=1.0,
                   xlab=bquote(~Log[2]~'fold change'~.(grp)~'/'~.(controlGrp)),
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
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
   "_volcanoPlot_chrX.pdf"), plot=p3,
   device="pdf",path=outPath,width=12,height=12,units="cm")
   } else {
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                             "_volcanoPlot_chrX.png"), plot=p3,
             device="png",path=outPath,width=12,height=12,units="cm")
   }


   idx<-resByChr$chr!="chrX"
   #pdf(file=paste0(outPath,"/plots/",fileNamePrefix, grp,
   #                "_volcanoPlot_autosomes.pdf"), width=8,height=6,paper="a4")
   sigUp<-sum(resByChr$padj[idx]<pThresh & resByChr$log2FoldChange[idx]>1)
   sigDown<-sum(resByChr$padj[idx]<pThresh & resByChr$log2FoldChange[idx]< -1)
   p4<-EnhancedVolcano(resByChr[idx,],
                   lab=rownames(resByChr[idx,]),
                   labSize=0.5,
                   labCol="#11111100",
                   x="log2FoldChange",y="padj",
                   selectLab=rownames(resByChr)[12366],
                   xlim=c(-5.5,5.5),
                   ylim=c(0,65),
                   title= paste0(grp," vs ",controlGrp,": autosomal genes"),
                   subtitle=NULL,
                   caption = paste0(nrow(resLFC), ' expressed genes. ',sigUp, " up, ",sigDown," down."),
                   captionLabSize = 12,
                   pCutoff=pThresh,
                   FCcutoff=1.0,
                   xlab=bquote(~Log[2]~'fold change'~.(grp)~'/'~.(controlGrp)),
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
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                             "_volcanoPlot_autosomes.pdf"), plot=p4,
             device="pdf",path=outPath,width=12,height=12,units="cm")
   } else {
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                      "_volcanoPlot_autosomes.png"), plot=p4,
          device="png",path=outPath,width=12,height=12,units="cm")
   }


   summaryByChr<-function(resLFC,pThresh,LFCthresh) {
      up<-resLFC[resLFC$padj < pThresh & resLFC$log2FoldChange > LFCthresh,]
      down<-resLFC[resLFC$padj < pThresh & resLFC$log2FoldChange < -LFCthresh, ]
      allChr<-as.data.frame(rbind(up=table(up$chr),down=table(down$chr)))
      allChr$autosomes<-rowSums(allChr[,1:5])
      allChr$total<-rowSums(allChr[,1:6])
      rownames(allChr)<-paste0(rownames(allChr),"_p",pThresh,"_lfc",LFCthresh)
      return(allChr)
   }


   sink(file=paste0(outPath,"/txt/", fileNamePrefix, grp,
                    "_logfile.txt"),append=TRUE, type="output")
   cat("Summary by Chr: \n")
   cat("\np=0.05, LFC=0: \n")
   print(summaryByChr(resLFC,pThres=0.05,LFCthresh=0))
   cat("\np=0.05, LFC=0.5: \n")
   print(summaryByChr(resLFC,pThres=0.05,LFCthresh=0.5))
   cat("\np=0.05, LFC=1: \n")
   print(summaryByChr(resLFC,pThres=0.05,LFCthresh=1))

   cat("\np=0.01, LFC=0: \n")
   print(summaryByChr(resLFC,pThres=0.01,LFCthresh=0))
   cat("\np=0.01, LFC=0.5: \n")
   print(summaryByChr(resLFC,pThres=0.01,LFCthresh=0.5))
   cat("\np=0.01, LFC=1: \n")
   print(summaryByChr(resLFC,pThres=0.01,LFCthresh=1))

   sink()
   #summary(resLFC,alpha=0.05)


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


   salmondc<-filterResults(results,padj=0.05,lfc=0.5,"gt","chrX", writeTable=F)
   salmondcgr<-metadata[metadata$wormbaseID %in% salmondc$wormbaseID]
   mcols(salmondcgr)<-cbind(mcols(salmondcgr),
                            salmondc[match(salmondcgr$wormbaseID,
                                           salmondc$wormbaseID),c(1:3)])
   salmondcgr
   lfcVal=0.5 # >log2(1.4) and <log2(1.5)
   padjVal=0.05
   saveRDS(salmondcgr,file=paste0(outPath,"/rds/",fileNamePrefix, grp,
                                  "_chrXup_lfc",
                                  formatC(lfcVal,format="e",digits=0),"_p",
                                  formatC(padjVal,format="e",digits=0),
                                  "_gr.rds"))
}



# Volcano - colour by other datasets --------------------------------------



oscillating<-read.delim(paste0(outPath,"/oscillatingGenes.tsv"),header=T,
                        stringsAsFactors=F)

hsUp<-readRDS("hsUp_garrigues2019.rds")
hsDown<-readRDS("hsDown_garrigues2019.rds")

padjVal=0.05
lfcVal=0.5

#grp=groupsOI[3]

for(grp in groupsOI){
   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

   #### oscillating genes
   keyvals<-rep('black', nrow(salmon))
   names(keyvals)<-rep('Other',nrow(salmon))
   idx<-salmon$wormbaseID %in% oscillating$WB_ID
   keyvals[idx]<-'red'
   names(keyvals)[idx]<-"Oscillating"
   sigUp<-sum(rowSums(cbind(salmon$padj< padjVal,
                            salmon$log2FoldChange>lfcVal,
                            keyvals=="red"), na.rm=T)==3, na.rm=T)
   sigDown<-sum(rowSums(cbind(salmon$padj< padjVal,
                              salmon$log2FoldChange<lfcVal,
                              keyvals=="red"), na.rm=T)==3, na.rm=T)
   p1<-EnhancedVolcano(salmon,
                       lab=salmon$publicID,
                       labSize=0.5,
                       labCol="#11111100",
                       x="log2FoldChange",
                       y="padj",
                       selectLab=salmon$publicID[12366],
                       xlim=c(-5.5,5.5),
                       ylim=c(0,65),
                       title= paste0(grp," vs ", controlGrp),
                       subtitle=NULL,
                       caption = paste0(sum(keyvals=="red"), ' oscillating genes (Meeuse 2020). ',sigUp, " up, ",sigDown," down."),
                       captionLabSize = 12,
                       pCutoff=padjVal,
                       FCcutoff=lfcVal,
                       xlab=bquote(~Log[2]~'fold change'~.(grp)~'/'~.(controlGrp)),
                       ylab=bquote(~-Log[10]~adjusted~italic(P)),
                       #.legend=c('NS','P & Log2 FC'),
                       #legendLabels=c('NS', expression(p-value<pThresh~and~log[2]~FC>1)),
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
                             "_volcanoPlot_oscillatingGenes.pdf"), plot=p1,
             device="pdf",path=outPath, width=12,height=12,units="cm")
   } else {
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                             "_volcanoPlot_oscillatingGenes.png"), plot=p1,
             device="png",path=outPath, width=12,height=12,units="cm")
   }



   #### heat shock genes
   myCols<-c("#11111100","red","red") # background, dataset1, dataset2
   keyvals<-rep(myCols[1], nrow(salmon))
   names(keyvals)<-rep('Other',nrow(salmon))
   idx<-salmon$wormbaseID %in% hsUp$WormBase.ID | salmon$wormbaseID %in% hsDown$WormBase.ID
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
                       title= paste0(grp," vs ", controlGrp),
                       subtitle=NULL,
                       caption = paste0(sum(keyvals==myCols[2]), ' heatshock genes (Garrigues 2019). ',sigUp, " up, ",sigDown," down."),
                       captionLabSize = 12,
                       pCutoff=padjVal,
                       FCcutoff=lfcVal,
                       xlab=bquote(~Log[2]~'fold change'~.(grp)~'/'~.(controlGrp)),
                       ylab=bquote(~-Log[10]~adjusted~italic(P)),
                       #.legend=c('NS','P & Log2 FC'),
                       #legendLabels=c('NS', expression(p-value<pThresh~and~log[2]~FC>1)),
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
                             "_volcanoPlot_hsGenes.pdf"), plot=p1,
             device="pdf",path=outPath, width=12,height=12,units="cm")
   } else {
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                             "_volcanoPlot_hsGenes.png"), plot=p1,
             device="png",path=outPath, width=12,height=12,units="cm")
   }
}


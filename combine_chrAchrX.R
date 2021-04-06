library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library(GenomicFeatures)
library(ggplot2)
library("RColorBrewer")
library(DESeq2)
library(tximport)
#library("PoiClaClu")
#library("pheatmap")
library(tidyr)
library(EnhancedVolcano)
#library(affy)
#library("gplots")
library(ggpubr)
library(plyr)
library(dplyr)


source("./variableSettings.R")
source("./functions.R")
metadata<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))
if(filterData){
  fileNamePrefix<-filterPrefix
}


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


fileList<-read.table(paste0(outPath,"/fastqList.txt"),stringsAsFactors=F,header=T)

sampleNames<-paste(fileList$sampleName, fileList$repeatNum, fileList$laneNum, sep="_")

fileNames<-paste0(outPath,"/salmon/mRNA/",sampleNames,"/quant.sf")

sampleTable<-data.frame(fileName=fileNames,sampleName=sampleNames,stringsAsFactors=F)

# extract the technical replicate variable
sampleTable$replicate=factor(fileList$repeatNum)
sampleTable$lane=factor(fileList$laneNum)

# extract the strain variable
sampleTable$strain<-factor(as.character(fileList$sampleName),levels=c("366","382","775","784"))
sampleTable$SMC<-sampleTable$strain
levels(sampleTable$SMC)<-c("wt","dpy26cs","kle2cs","scc1cs")

controlGrp<-levels(sampleTable$SMC)[1] # control group
groupsOI<-levels(sampleTable$SMC)[-1] # groups of interest to contrast to control

fileNamePrefix<-filterPrefix



# # Create metadata object --------------------------------------------------
# ###############################################################-
# ### create metadata
# ###############################################################-
#
# if(!file.exists(paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.",
#                        genomeVer, ".annotations.sqlite"))){
#   dir.create(paste0(genomeDir,"/annotations"),recursive=T)
#   system(paste0("curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/gff/c_elegans.PRJNA13758.",
#                 genomeVer,".annotations.gff3.gz -o ",genomeDir,
#                 "/annotations/c_elegans.PRJNA13758.",genomeVer,
#                 ".annotations.gff3.gz"))
#
#   system(paste0("gunzip ",genomeDir,"/annotations/c_elegans.PRJNA13758.",
#                 genomeVer,".annotations.gff3.gz"))
#   si<-seqinfo(Celegans)
#   genome(si)<-genomeVer
#   seqnames(si)<-gsub("M","MtDNA",gsub("chr","",seqnames(si)))
#   wstxdb<-makeTxDbFromGFF(file=paste0(genomeDir,
#                                       "/annotations/c_elegans.PRJNA13758.",
#                                       genomeVer,".annotations.gff3"),
#                           format="gff3",organism="Caenorhabditis elegans",
#                           chrominfo=si)
#
#   saveDb(wstxdb,paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.", genomeVer,
#                        ".annotations.sqlite"))
#   file.remove(paste0(genomeDir,"/annotations/c_elegans.PRJNA13758.",
#                      genomeVer, ".annotations.gff3"))
# }
#
# # load a txdb of wormbase data and create a tx2gene object
# txdb<-loadDb(paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.", genomeVer,
#                     ".annotations.sqlite"))
# k <- keys(txdb, keytype = "TXNAME")
# tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
#
# columns(txdb)
# keytypes(txdb)
# TxptByGene<-transcriptsBy(txdb, by = "gene")
# length(TxptByGene)
#
# geneGR<-unlist(range(TxptByGene))
# mcols(geneGR)$wormbase<-names(geneGR)
# genedf<-as.data.frame(geneGR)
#
#
# # download gene id data from simplemine: https://wormbase.org/tools/mine/simplemine.cgi
# # for entrez ids, load wormbaseID column into https://david.ncifcrf.gov/conversion.jsp
# geneIDs<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/simplemine_WS278_geneID.txt")
# david<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/david_wbid2entrez_WS278.txt")
#
# metadata<-inner_join(geneIDs, genedf,by=c("WormBase.Gene.ID"="wormbase")) %>%
#   dplyr::select(WormBase.Gene.ID,Public.Name,Sequence.Name,seqnames,start, end, strand) %>%
#   collect %>% GenomicRanges::GRanges()
#
# names(mcols(metadata))<-c("wormbaseID","publicID","sequenceID")
#
# i<-which(metadata$wormbaseID %in% david$From)
# j<-match(metadata$wormbaseID[i],david$From)
# metadata$entrezID<-NA
# metadata$entrezID[i]<-david$To[j]
#
# #seqinfo(metadata)<-wbseqinfo
# seqlevelsStyle(metadata)<-"ucsc"
# seqinfo(metadata)<-ce11seqinfo
# metadata<-sort(metadata)
#
# saveRDS(metadata,paste0(outPath,"/wbGeneGR_WS275.rds"))
#
# ###############################################################-
# # Import into DESeq2 ------------------------------------------------------
# ###############################################################-
#
# # import the count matrices
# txi<-tximport(sampleTable$fileName,type="salmon",tx2gene=tx2gene)
#
# # read samples into DESeq2
# dds <- DESeqDataSetFromTximport(txi=txi,
#                                 colData=sampleTable,
#                                 design=~replicate+lane+SMC)
#
#
#
# ###############################################################-
# ### DESeq2 differential expression analysis (using negative binomial distribution)
# ###############################################################-
#
# #dds<-collapseReplicates(dds,groupby=samples$sampleID,renameCols=T)
# #dds <- DESeq(dds)
# # This function performs a default analysis through the steps:
# #   1. estimation of size factors: estimateSizeFactors
# #   2. estimation of dispersion: estimateDispersions
# #   3. Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
# #   returns a DESeqDataSet object
#
# idx<-match(rownames(dds),metadata$wormbaseID)
# # add gene and chormosome names as metadata
# featureData <- data.frame(gene=rownames(dds),
#                           chr=as.vector(seqnames(metadata))[idx]) #,
#
# rowData(dds) <- DataFrame(mcols(dds), featureData)
#
# #remove unmapped or mitochondrial genes
# idx<-is.na(rowData(dds)$chr) | rowData(dds)$chr %in% c("MtDNA","chrM")
# dds<-dds[!idx,]
#
# #prefilter rows with less than 10 reads in total
# dds<-dds[rowSums(counts(dds)) >= 10,]
# print(paste0(dim(dds)[1], " genes with >=10 reads total"))
#
# ###################-
# ####### filter genes-------------------------
# #######-###########-
# if(filterData){
#   # remove filtered genes
#   idx<-rowData(dds)$gene %in% toFilter
#   dds<-dds[!idx,]
#
#   fileNamePrefix=filterPrefix
# }
#
# # do not run statistical testing
# dds<-DESeq(dds)
#
# saveRDS(dds,file=paste0(outPath,"/rds/dds_object.rds"))




######################-
## combine chrX and chrA #####
######################-

chrXprefix<-paste0("/../SMC_RNAseq_prefiltCyc2x/rds/p",padjVal,"_lfc",lfcVal,"/preFiltOsc2x_")
chrAprefix<-paste0("/../SMC_RNAseq_prefiltCyc2xChrA/rds/p",padjVal,"_lfc",lfcVal,"/preFiltOsc2xChrA_")

for (grp in groupsOI){
  #grp=groupsOI[1]
  Xdata<-readRDS(paste0(outPath, chrXprefix, grp, "_DESeq2_fullResults.rds"))
  table(Xdata$chr)
  Adata<-readRDS(paste0(outPath, chrAprefix, grp, "_DESeq2_fullResults.rds"))
  table(Adata$chr)
  resLFC<-rbind(Adata,Xdata[Xdata$chr=="chrX",])
  saveRDS(resLFC,paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))
  #export csv with ordered results
  write.csv(resLFC[order(resLFC$padj),],
            file=paste0(outPath,"/txt/", fileNamePrefix,grp,
                        "_DESeq2_resultsTable.csv"),
            quote=F,row.names=F)

  # remove NAs from chr (unmapped or mtDNA) and padj (below filter threshold) columns
  idx<-is.na(resLFC$chr) | is.na(resLFC$padj)
  #res<-res[!idx,]
  resLFC<-resLFC[!idx,]


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
  score<-oldf %>% group_by(queryHits) %>% dplyr::summarise(score=mean(scorePerQuery))
  forBW$score<-score$score
  export(forBW,paste0(outPath,"/tracks/",fileNamePrefix,grp,
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
  export(forBed,paste0(outPath,"/tracks/",fileNamePrefix,grp,
                       "_wt_lfc_p",gsub("^0.","",padjVal),".bedGraph"),
         format="bedGraph")


  export(forBed,paste0(outPath,"/tracks/",fileNamePrefix,grp,
                       "_wt_lfc_p",gsub("^0.","",padjVal),".bed"),
         format="bed")





  # MAplots -----------------------------------------------------------------
  ###########-
  # MAplot ALL genes
  ############-

  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,grp,
                  "_MAplots_results.pdf"), width=5,height=5,paper="a4")

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
                                             fileNamePrefix, grp,
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

  sink(file=paste0(outPath,"/txt/",fileNamePrefix,grp,
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
  pdf(file=paste0(outPath,"/plots/",fileNamePrefix, grp,
                  "_boxPlots_expnByChrType.pdf"), width=5,height=5,paper="a4")

  idx<-resLFC$log2FoldChange!=0
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
                      title= paste0(grp," vs ", controlGrp),
                      subtitle=NULL,
                      caption = paste0(nrow(resLFC), ' expressed genes. ',sigUp,
                                       " up, ",sigDown," down."),
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
                      title= paste0(grp," vs ",controlGrp),
                      subtitle=NULL,
                      caption = paste0(nrow(resLFC), ' expressed genes. ',sigUp, " up, ",sigDown," down."),
                      captionLabSize = 12,
                      pCutoff=padjVal,
                      FCcutoff=lfcVal,
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

  if(length(chrXgenes)>0) {
    idx<-resByChr$chr=="chrX"
    #pdf(file=paste0(outPath,"/plots/",fileNamePrefix, grp,
    #                "_volcanoPlot_chrX.pdf"), width=8,height=6,paper="a4")
    sigUp<-sum(resByChr$padj[idx]<padjVal & resByChr$log2FoldChange[idx]>lfcVal)
    sigDown<-sum(resByChr$padj[idx]<padjVal & resByChr$log2FoldChange[idx]< -lfcVal)
    p3<-EnhancedVolcano(resByChr[idx,],
                        lab=rownames(resByChr[idx,]),
                        x="log2FoldChange",y="padj",
                        selectLab=rownames(resByChr)[12366],
                        xlim=c(-5.5,5.5),
                        ylim=c(0,65),
                        title= paste0(grp," vs ",controlGrp,": chrX genes"),
                        subtitle=NULL,
                        caption = paste0(sum(resByChr$chr=="chrX"), ' expressed genes. ',sigUp, " up, ",sigDown," down."),
                        captionLabSize = 12,
                        pCutoff=padjVal,
                        FCcutoff=lfcVal,
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
  }

  idx<-resByChr$chr!="chrX"
  #pdf(file=paste0(outPath,"/plots/",fileNamePrefix, grp,
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
                      title= paste0(grp," vs ",controlGrp,": autosomal genes"),
                      subtitle=NULL,
                      caption = paste0(sum(resByChr$chr!="chrX"), ' expressed genes. ',sigUp, " up, ",sigDown," down."),
                      captionLabSize = 12,
                      pCutoff=padjVal,
                      FCcutoff=lfcVal,
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


  summaryByChr<-function(resLFC,padj,lfc) {
    up<-resLFC[resLFC$padj < padjVal & resLFC$log2FoldChange > lfcVal,]
    down<-resLFC[resLFC$padj < padjVal & resLFC$log2FoldChange < -lfcVal, ]
    allChr<-as.data.frame(rbind(up=table(up$chr),down=table(down$chr)))
    allChr$autosomes<-rowSums(allChr[,1:5])
    allChr$total<-rowSums(allChr[,1:6])
    rownames(allChr)<-paste0(rownames(allChr),"_p",padjVal,"_lfc",lfcVal)
    return(allChr)
  }


  sink(file=paste0(outPath,"/txt/", fileNamePrefix, grp,
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
    salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))
    salmondc<-filterResults(salmon,padj=padjVal,lfc=lfcVal,"gt","chrX", writeTable=F)
    salmondcgr<-metadata[metadata$wormbaseID %in% salmondc$wormbaseID]
    mcols(salmondcgr)<-cbind(mcols(salmondcgr),
                             salmondc[match(salmondcgr$wormbaseID,
                                            salmondc$wormbaseID),c(1:3)])
    salmondcgr
    saveRDS(salmondcgr,file=paste0(outPath,"/rds/",fileNamePrefix, grp,
                                   "_chrXup_lfc", lfcVal,"_p",
                                   padjVal, "_gr.rds"))
  }
}



# Volcano - colour by other datasets --------------------------------------

oscillating<-read.delim(paste0(outPath,"/publicData/oscillatingGenes.tsv"),header=T,
                        stringsAsFactors=F)
latorre<-read.delim(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv")) #3235
hsUp<-readRDS(paste0(outPath,"/publicData/hsUp_garrigues2019.rds"))
hsDown<-readRDS(paste0(outPath,"/publicData/hsDown_garrigues2019.rds"))

amplicons<-readRDS(paste0(outPath,"/otherData/ampliconMaxTSSgr.RDS"))


#grp=groupsOI[3]

for(grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

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
                      title= paste0(grp," vs ", controlGrp),
                      subtitle=NULL,
                      caption = paste0(sum(keyvals!=bkgrnd), ' oscillating genes. ',sigUp, " up, ",sigDown," down."),
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
                           "_volcanoPlot_oscillatingGenes.pdf"), plot=p1,
           device="pdf",path=outPath, width=12,height=12,units="cm")
  } else {
    ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                           "_volcanoPlot_oscillatingGenes.png"), plot=p1,
           device="png",path=outPath, width=12,height=12,units="cm")
  }


  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))



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
                      title= paste0(grp," vs ", controlGrp),
                      subtitle=NULL,
                      caption = paste0(sum(keyvals==myCols[2]), ' heatshock genes (Garrigues 2019). ',sigUp, " up, ",sigDown," down."),
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
                           "_volcanoPlot_hsGenes.pdf"), plot=p1,
           device="pdf",path=outPath, width=12,height=12,units="cm")
  } else {
    ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
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
                      title= paste0(grp," vs ", controlGrp),
                      subtitle=NULL,
                      caption = paste0(sum(keyvals!="black"), ' amplicon genes: ',sigUp, " up, ",sigDown," down."),
                      captionLabSize = 12,
                      pCutoff=padjVal,
                      FCcutoff=lfcVal,
                      xlab=bquote(~Log[2]~'fold change'~.(grp)~'/'~.(controlGrp)),
                      ylab=bquote(~-Log[10]~adjusted~italic(P)),
                      legendPosition = 'top',
                      legendLabSize = 12,
                      legendIconSize = 3.0,
                      axisLabSize=14,
                      colCustom=keyvals,
                      colAlpha=0.5,
                      pointSize = 1.0)

  if(plotPDFs==T){
    ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                           "_volcanoPlot_ampliconGenes.pdf"), plot=p1,
           device="pdf",path=outPath, width=12,height=12,units="cm")
  } else {
    ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                           "_volcanoPlot_ampliconGenes.png"), plot=p1,
           device="png",path=outPath, width=12,height=12,units="cm")
  }

}



########################-
## ECDF of data -----
########################-

sigTables<-list()
localPadj=0.05
localLFC=0
for (grp in groupsOI){
  print(grp)
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))
  print(dim(salmon))
  print(sum(is.na(salmon$log2FoldChange)))
  #salmon$expressed<-sum(salmon$baseMean>10)
  sigTables[[prettyGeneName(grp)]]<-as.data.frame(salmon[salmon$baseMean>10,])
  print(dim(sigTables[[prettyGeneName(grp)]]))
}

# check if datasets have chrX genes included
includeChrX<-"chrX" %in% unlist(lapply(sigTables,"[","chr"))

SMC<-rep(names(sigTables),lapply(sigTables,nrow))
sig<-do.call(rbind,sigTables)
sig$SMC<-SMC
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


dd<-plyr::ddply(sig,.(SMC,upVdown,XvA),transform,
                ecd=ecdf(abs(log2FoldChange))(abs(log2FoldChange)))

dd1<-sig %>% dplyr::group_by(SMC,upVdown,XvA) %>%
  mutate(ecd=ecdf(abs(log2FoldChange))(abs(log2FoldChange)))

sig %>% dplyr::group_by(SMC,XvA) %>%
  mutate(countOnChr=n()) %>% group_by(SMC,XvA,upVdown) %>%
  mutate(countOnChr=unique(countOnChr),countInGrp=n(),
         fractionInGrp=countInGrp/countOnChr) %>%
  summarise(ecd=1-ecdf(abs(log2FoldChange))(c(0)),
            fractionInGrp=unique(fractionInGrp),
            countOnChr=unique(countOnChr),
            countInGrp=unique(countInGrp),
            countPosInGrp=ecd*countInGrp,
            percentPosInGrp=ecd*countInGrp/countOnChr)

sig %>% dplyr::group_by(SMC,XvA) %>%
  mutate(countOnChr=n()) %>% group_by(SMC,XvA,upVdown) %>%
  mutate(countOnChr=unique(countOnChr),countInGrp=n(),
         fractionInGrp=countInGrp/countOnChr) %>%
  summarise(ecd=1-ecdf(abs(log2FoldChange))(c(0)),
            fractionInGrp=unique(fractionInGrp),
            countOnChr=unique(countOnChr),
            countInGrp=unique(countInGrp),
            countPosInGrp=ecd*countInGrp,
            percentPosInGrp=ecd*countInGrp/countOnChr,
            sig=sum(abs(log2FoldChange)>localLFC & padj<localPadj))

ss<- sig %>% filter(abs(log2FoldChange)>localLFC,padj<localPadj) %>%
  group_by(SMC,XvA,upVdown) %>% summarise(count=n())

 ss


p<-ggplot(dd1, aes(x=abs(log2FoldChange),y=ecd,color=SMC,linetype=XvA)) +
  geom_line(size=1)+ facet_wrap(vars(upVdown),nrow=2)+
  theme_classic() + xlim(c(0,1.5)) +
  xlab("Absolute log2 fold change")+ylab("Fraction genes rejected")
p

#stat_ecdf(aes(colour=SMC,linetype=XvA),alpha=0.7)
p1<-p+geom_vline(aes(xintercept = 0.5), color="grey") +
  annotate("text",label="0.5",size=3, x=0.5, y=0,hjust=-0.05,color="grey") +
  geom_vline(aes(xintercept = 0.25), color="grey") +
  annotate("text",label="0.25",size=3, x=0.25, y=0,hjust=-0.05,color="grey")
p1

if(plotPDFs==T){
  ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix,
                         "lfcValueCDF.pdf"), plot=p1,
         device="pdf",path=outPath, width=10,height=10,units="cm")
} else {
  ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                         "lfcValueCDF.png"), plot=p1,
         device="png",path=outPath, width=10,height=10,units="cm")
}

dd %>% group_by(SMC,XvA) %>% mutate(expressed=n()) %>% group_by(SMC,XvA,upVdown) %>%
  summarise(qnt25=1-ecdf(abs(log2FoldChange))(0.25),
            qnt50=1-ecdf(abs(log2FoldChange))(0.5),
            count=n(), expressed=unique(expressed)) %>%
  mutate(percent0.25=100*count*qnt25/expressed,
         percent0.5=100*count*qnt50/expressed)


table(sig$SMC,sig$XvA,sig$upVdown)

library(rtracklayer)
#library(GenomicInteractions)
library(ggplot2)
#library(EnhancedVolcano)
library(BSgenome.Celegans.UCSC.ce11)
library(zoo)
library(dplyr)
library(ggpubr)
library(genomation)
library(seqplots)
library(RColorBrewer)
library(ggpubr)

source("functions.R")
source("./variableSettings.R")

scriptName <- "compareHICfeatures"
print(scriptName)

if(filterData){
  fileNamePrefix<-filterPrefix
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
} else {
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
}

makeDirs(outPath,dirNameList=paste0(c("plots/","tracks/"),
                                    paste0(dirname(fileNamePrefix),"/",
                                           scriptName)))


####### Supplementary functions---------------

#' Collect counts of significantly changed genes per chromosome and per compartment
#'
#' @param listgr List of GRanges for different RNAseq with PCA data in mcols
#' @param namePCAcol Name of PCA column in the listgr object
#' @param padjVal adjusted p value to use as threshold
#' @param lfcVal Log2 fold change value to use as threshold
#' @result Table of significant genes with significant up and down regulated genes
#' counted by chromosome for each sample
#' @export
processCountsPerChr<-function(listgr,namePCAcol,padjVal=0.05,lfcVal=0.5){
  # genes that change significantly
  sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                  padj=padjVal,lfc=lfcVal,direction="both")

  # count genes by category (chr & A/B)
  dfl<-lapply(sigList, function(x){x%>% dplyr::group_by(seqnames, get(namePCAcol)) %>% tally()})
  bgCount<-lapply(lapply(listgr,as.data.frame), function(x){x%>% dplyr::group_by(seqnames,get(namePCAcol)) %>% tally()})
  # add name of SMC protein
  dfl<-do.call(rbind, mapply(cbind,dfl,"SMC"=names(dfl),SIMPLIFY=F))
  bgCount<-do.call(rbind, mapply(cbind,bgCount,"SMC"=names(bgCount),SIMPLIFY=F))
  names(dfl)<-c("seqnames",namePCAcol,"n","SMC")
  names(bgCount)<-c("seqnames",namePCAcol,"n","SMC")
  dfl$seqnames<-gsub("chr","",dfl$seqnames)
  bgCount$seqnames<-gsub("chr","",bgCount$seqnames)

  # do left join to make sure dfl has all the categories required
  dfl<-left_join(bgCount,dfl,by=c("seqnames",namePCAcol,"SMC"),suffix=c("_total",""))
  dfl$n[is.na(dfl$n)]<-0
  dfl$Frac<-dfl$n/dfl$n_total
  return(dfl)
}

#' Plot counts of significantly changed genes per chromosome and per compartment
#'
#' @param df Data frame with counts of significant genes per chromosome and compartment
#' @param namePCAcol Name of PCA column in the df object
#' @param padjVal Name of eigen vector to use in plot title
#' @result Plot of significant genes with significant changed genes
#' counted by chromosome for each sample
#' @export
plotCountsPerChrPerCompartment<-function(df,namePCAcol,namePCA){
  ymax=max(df$n)
  p<-ggplot(df,aes(x=seqnames,y=n,group=get(namePCAcol))) +
    geom_bar(stat="identity", position=position_dodge(), aes(fill=get(namePCAcol))) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
    xlab("chr")+ylab("Number of genes") +
    theme(legend.title = element_text(size=10)) +
    ggtitle(paste0("Significantly changed genes per chromosome by ",namePCA," compartment")) + labs(fill=namePCA)
  return(p)
}


#' Plot fractions of significantly changed genes per chromosome and per compartment
#'
#' @param df Data frame with fractions of significant genes per chromosome and compartment
#' @param namePCAcol Name of PCA column in the df object
#' @param padjVal Name of eigen vector to use in plot title
#' @result Plot of significant genes with significant changed genes
#' counted by chromosome for each sample presented as fraction
#' @export
plotFractionPerChrPerCompartment<-function(df,namePCAcol,namePCA){
  p<-ggplot(df,aes(x=seqnames,y=Frac,group=get(namePCAcol))) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=get(namePCAcol))) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
    xlab("chr")+ylab("Number of genes") +
    ggtitle(paste0("Fraction changed genes per chromosome by ",namePCA," compartment")) + labs(fill=namePCA)
  return(p)
}

#' Collect counts of significantly up/down genes per chromosome and per compartment
#'
#' @param listgr List of GRanges for different RNAseq with PCA data in mcols
#' @param namePCAcol Name of PCA column in the listgr object
#' @param padjVal adjusted p value to use as threshold
#' @param lfcVal Log2 fold change value to use as threshold
#' @result Table of significant genes with significant up and down regulated genes
#' counted by chromosome for each sample
#' @export
processUpDownCountsPerChrPerCompartment<-function(listgr, namePCAcol, namePCA, padjVal=0.05, lfcVal=0.5){
  # upregulated genes
  sigListUp<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                    padj=padjVal,lfc=lfcVal,direction="gt")
  # count genes by category (chr & A/B)
  dflUp<-lapply(sigListUp, function(x){x%>% dplyr::group_by(seqnames,get(namePCAcol)) %>% tally()})
  # add name of SMC protein
  dflUp<-do.call(rbind, mapply(cbind,dflUp,"SMC"=names(dflUp),SIMPLIFY=F))
  names(dflUp)<-c("seqnames",namePCAcol,"n","SMC")

  #get bgCount for full category list
  bgCount<-lapply(lapply(listgr,as.data.frame), function(x){x%>% dplyr::group_by(seqnames,get(namePCAcol)) %>% tally()})
  # add name of SMC protein
  bgCount<-do.call(rbind, mapply(cbind,bgCount,"SMC"=names(bgCount),SIMPLIFY=F))
  names(bgCount)<-c("seqnames",namePCAcol,"n","SMC")

  # do left join to make sure dfl has all the categories required
  dflUp<-dplyr::left_join(bgCount,dflUp,by=c("seqnames",namePCAcol,"SMC"),suffix=c("_total",""))
  dflUp$n[is.na(dflUp$n)]<-0
  dflUp$expression<-"up"

  # downregulated genes
  sigListDown<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                      padj=padjVal,lfc= -lfcVal, direction="lt")
  dflDown<-lapply(sigListDown, function(x){x%>% dplyr::group_by(seqnames,get(namePCAcol)) %>% tally()})
  # add name of SMC protein
  dflDown<-do.call(rbind, mapply(cbind,dflDown,"SMC"=names(dflDown),SIMPLIFY=F))
  names(dflDown)<-c("seqnames",namePCAcol,"n","SMC")

  # do left join to make sure dfl has all the categories required
  dflDown<-left_join(bgCount,dflDown,by=c("seqnames",namePCAcol,"SMC"),suffix=c("_total",""))
  dflDown$n[is.na(dflDown$n)]<-0
  dflDown$expression<-"down"

  dfl<-rbind(dflUp,dflDown)
  dfl$expression<-factor(dfl$expression,levels=c("up","down"))
  dfl$seqnames<-gsub("chr","",dfl$seqnames)
  return(dfl)
}

#' Plot counts of significantly up/down genes per compartment
#'
#' @param df Data frame with fractions of significant genes per chromosome and compartment
#' @param namePCAcol Name of PCA column in the df object
#' @param namePCA Name of eigen vector to use in plot title
#' @param compartment Name of compartment to plot ("A" or "B")
#' @result Plot of significant genes with significant changed genes
#' counted by chromosome for each sample presented as fraction
#' @export
plotUpDownByCompartment<-function(df,namePCAcol,namePCA,compartment){
  yminmax=c(0,max(df$n))
  p<-ggplot(df[df[,namePCAcol]==compartment,],aes(x=seqnames,y=n,group=expression)) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=expression)) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
    xlab("chr")+ylab("Number of genes") +
    ggtitle(paste0("Up/down regulated in ",namePCA," ",compartment,
                   " compartment"))
  return(p)
}

#' Collect LFC of significantly changed genes per chromosome and per compartment
#'
#' @param listgr List of GRanges for different RNAseq with PCA data in mcols
#' @param namePCAcol Name of PCA column in the listgr object
#' @param padjVal adjusted p value to use as threshold
#' @param lfcVal Log2 fold change value to use as threshold
#' @result Table of significant genes with LFC of significant up and down regulated genes
#' for each sample
#' @export
processLFCbyCompartment<-function(listgr,namePCAcol,padjVal=0.05,lfcVal=0){
  # upregulated
  sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                  padj=padjVal,lfc=lfcVal,direction="gt")

  sigList<-lapply(sigList, "[", ,c(namePCAcol,"log2FoldChange"))
  for(g in names(sigList)){ sigList[[g]]$SMC<-g }
  sigList<-do.call(rbind,sigList)
  sigList<-sigList[!is.na(sigList[,namePCAcol]),]
  sigTbl<-sigList
  sigTbl$updown<-"up"

  # downregulated
  sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                  padj=padjVal, lfc= -lfcVal, direction="lt")

  sigList<-lapply(sigList, "[", ,c(namePCAcol,"log2FoldChange"))
  for(g in names(sigList)){ sigList[[g]]$SMC<-g }
  sigList<-do.call(rbind,sigList)
  sigList<-sigList[!is.na(sigList[,namePCAcol]),]
  sigList$updown<-"down"
  sigTbl<-rbind(sigTbl,sigList)
  sigTbl[,namePCAcol]<-as.factor(sigTbl[,namePCAcol])
  sigTbl$updown<-factor(sigTbl$updown, levels=c("up","down"))
  return(sigTbl)
}

give.n <- function(x){
  return(c(y = 0, label = length(x)))
  # experiment with the multiplier to find the perfect position
}


#' Plot LFC of significantly up/down genes per compartment
#'
#' @param df Data frame with LFC of significant genes in A/B compartments
#' @param namePCAcol Name of PCA column in the df object
#' @param namePCA Name of eigen vector to use in plot title
#' @param compartment Name of compartment to plot ("A" or "B")
#' @result Plot of significant genes with significant changed genes
#' counted by chromosome for each sample presented as fraction
#' @export
plotLFCbyCompartment<-function(df,namePCAcol,namePCA){
  #df<-sigTbl[sigTbl[,namePCAcol]==compartment,]
  #df<-sigTbl
  yminmax=c(0,median(abs(df$log2FoldChange))+quantile(abs(df$log2FoldChange))[4]*2)
  p<-ggplot(df,aes(x=get(namePCAcol),y=abs(log2FoldChange),col=updown,fill=updown)) +
    geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2),
                 outlier.shape=NA) +
    facet_grid(cols=vars(SMC)) +
    ggtitle(paste0("Significantly changed genes by ",namePCA," compartment")) +
    theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
    scale_y_continuous(limits = yminmax) +
    scale_color_grey(start=0.2,end=0.2,guide="none") +
    stat_summary(fun.data = give.n, geom = "text",size=3,angle=90,
                 position=position_dodge(width=0.6))
  return(p)
}

#############################-
##### subsitute for getREF function from seqplots thaht has an unfixed bug.
##### Fix comes form  https://github.com/Przemol/seqplots/issues/58
#' Get reference genome
#'
#' @param genome The filename of FASTA file or genome code for BSgenome
#'
#' @return \code{DNAStringSet}
#'
#' @export
#'
getREF <- function(genome) {
  if( file.exists(file.path(Sys.getenv('root'), 'genomes', genome)) ) {
    REF <- Biostrings::readDNAStringSet( file.path(Sys.getenv('root'), 'genomes', genome) )
    names(REF) <- gsub(' .+', '', names(REF))
  } else {

    GENOMES <- BSgenome::installed.genomes(
      splitNameParts=TRUE)$genome
    if( length(GENOMES) )
      names(GENOMES) <- gsub('^BSgenome.', '', BSgenome::installed.genomes())
    if( !length(GENOMES) ) stop('No genomes installed!')

    pkg <- paste0('BSgenome.', names(GENOMES[GENOMES %in% genome]))[[1]]
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE, quietly=TRUE)
    )
    REF <- get(pkg)
  }
  return(REF)
}

assignInNamespace("getREF",getREF,ns="seqplots")

# why is pca1 just splitting the chromosome in two?
# mm<-matrix(data=rep(0,100),nrow=10)
# for (i in 1:10) {mm[i,i]<-10+rnorm(1)}
# for (i in 2:10) {mm[i,i-1]<-5+rnorm(1)}
# for (i in 2:10) {mm[i-1,i]<-5+rnorm(1)}
# for (i in 1:9) {mm[i,i+1]<-5+rnorm(1)}
# for (i in 1:9) {mm[i+1,i]<-5+rnorm(1)}
#
# pcaRes<-prcomp(mm)
# pcaRes$x[,1]
#
# plot(1:10,pcaRes$x[,1])
# lines(1:10,pcaRes$x[,1])
# abline(h=0)


# AB compartments - N2 ----------------------------------------------------

####
## AB compartments
####

####
## N2 compartments
####
pca1<-import.bw(paste0(outPath,"/otherData/N2_merge_2000.oriented_E1.vecs.bw"))
pca2<-import.bw(paste0(outPath,"/otherData/N2_merge_2000.oriented_E2.vecs.bw"))

listgr<-NULL
for (grp in useContrasts){
  #grp=useContrasts[1]
  salmon<-readRDS(file=paste0(paste0(outPath,"/rds/",fileNamePrefix,
                  contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))

  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

  salmongr<-sort(salmongr)
  salmongr<-assignGRtoAB(salmongr,pca1,grName=grp,pcaName="E1")
  salmongr<-assignGRtoAB(salmongr,pca2,grName=grp,pcaName="E2")
  listgr[[grp]]<-salmongr
}



####-
## AB compartment by chromosome-----
####-
######### N2 first Eigenvector ------


pcaSource="N2"
dfl<-processCountsPerChr(listgr,namePCAcol="E1_compartment")

p1<-plotCountsPerChrPerCompartment(dfl, namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"))
p1a<-plotFractionPerChrPerCompartment(dfl, namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"))

dfl<-processUpDownCountsPerChrPerCompartment(listgr, namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"))

p2<-plotUpDownByCompartment(dfl,namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"), compartment="A")
p3<-plotUpDownByCompartment(dfl,namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"), compartment="B")


p<-ggpubr::ggarrange(p1,p1a,p2,p3,ncol=1,nrow=2)
ggpubr::ggexport(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                                   "ABcomp_",pcaSource,"-E1_countsPerChr_padj",
                                   padjVal,"_lfc", lfcVal,".pdf"),
                 device="pdf",width=10,height=5, units="cm")

########## N2 second eigenvector ------

dfl<-processCountsPerChr(listgr,namePCAcol="E2_compartment")

p1<-plotCountsPerChrPerCompartment(dfl, namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"))
p1a<-plotFractionPerChrPerCompartment(dfl, namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"))

dfl<-processUpDownCountsPerChrPerCompartment(listgr, namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"))

p2<-plotUpDownByCompartment(dfl,namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"), compartment="A")
p3<-plotUpDownByCompartment(dfl,namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"), compartment="B")

p<-ggpubr::ggarrange(p1,p1a,p2,p3,ncol=1,nrow=2)
ggpubr::ggexport(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                                    "ABcomp_",pcaSource,"-E2_countsPerChr_padj",
                                    padjVal,"_lfc", lfcVal,".pdf"),
                  device="pdf",width=10,height=5, units="cm")


####-
## AB comp LFC-----
####-
#### first N2 eigenvector-----

pcaSource="N2"
localLFC=0
sigTbl<-processLFCbyCompartment(listgr,namePCAcol="E1_compartment", padjVal=padjVal,
                        lfcVal=localLFC)

p1<-plotLFCbyCompartment(sigTbl,namePCAcol="E1_compartment", namePCA=paste(pcaSource,"E1"))

ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                "ABcomp_",pcaSource,"-E1_LFC_padj",
                          padjVal,"_lfc", localLFC,".pdf"),
                plot=p1, device="pdf",width=29,height=19,units="cm")

# test significance of LFC
for(s in unique(sigTbl$SMC)){
  print(s)
  sigTblss<-sigTbl[sigTbl$SMC==s,]
  resAOV<-aov(abs(log2FoldChange)~updown*E1_compartment,data=sigTblss)
  print(summary(resAOV))
  #  TukeyHSD(resAOV)
}


### second N2 eigenvector------

localLFC=0
sigTbl<-processLFCbyCompartment(listgr,namePCAcol="E2_compartment", padjVal=padjVal,
                                lfcVal=localLFC)

p1<-plotLFCbyCompartment(sigTbl,namePCAcol="E2_compartment", namePCA=paste(pcaSource,"E2"))

ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                "ABcomp_",pcaSource,"-E2_LFC_padj",
                                padjVal,"_lfc", localLFC,".pdf"),
                plot=p1, device="pdf",width=29,height=19,units="cm")




####-
## AB compartment by chromosome 366-----
####-

####
## 366 compartments
####
pca1<-import.bw(paste0(outPath,"/otherData/366_merge_2000.oriented_E1.vecs.bw"))
pca2<-import.bw(paste0(outPath,"/otherData/366_merge_2000.oriented_E2.vecs.bw"))


listgr<-NULL
for (grp in useContrasts){
  #grp=useContrasts[1]
  salmon<-readRDS(file=paste0(paste0(outPath,"/rds/",fileNamePrefix,
                                     contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))

  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

  salmongr<-sort(salmongr)
  salmongr<-assignGRtoAB(salmongr,pca1,grName=grp,pcaName="E1")
  salmongr<-assignGRtoAB(salmongr,pca2,grName=grp,pcaName="E2")
  listgr[[grp]]<-salmongr
}

pcaSource="366"

### 366 first Eigenvector ------
dfl<-processCountsPerChr(listgr,namePCAcol="E1_compartment")

p1<-plotCountsPerChrPerCompartment(dfl, namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"))
p1a<-plotFractionPerChrPerCompartment(dfl, namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"))

dfl<-processUpDownCountsPerChrPerCompartment(listgr, namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"))

p2<-plotUpDownByCompartment(dfl,namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"), compartment="A")
p3<-plotUpDownByCompartment(dfl,namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"), compartment="B")


p<-ggpubr::ggarrange(p1,p1a,p2,p3,ncol=1,nrow=2)
ggpubr::ggexport(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                                   "ABcomp_", pcaSource, "-E1_countsPerChr_padj",
                                   padjVal,"_lfc", lfcVal,".pdf"),
                 device="pdf",width=10,height=5, units="cm")


### 366 second eigenvector ------
dfl<-processCountsPerChr(listgr,namePCAcol="E2_compartment")

p1<-plotCountsPerChrPerCompartment(dfl, namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"))
p1a<-plotFractionPerChrPerCompartment(dfl, namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"))

dfl<-processUpDownCountsPerChrPerCompartment(listgr, namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"))

p2<-plotUpDownByCompartment(dfl,namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"), compartment="A")
p3<-plotUpDownByCompartment(dfl,namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"), compartment="B")

p<-ggpubr::ggarrange(p1,p1a,p2,p3,ncol=1,nrow=2)
ggpubr::ggexport(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                                   "ABcomp_",pcaSource,"-E2_countsPerChr_padj",
                                   padjVal,"_lfc", lfcVal,".pdf"),
                 device="pdf",width=10,height=5, units="cm")



####-
## AB comp LFC-----
####-

localLFC=0
pcaSource="366"
### first 366 eigenvector------
sigTbl<-processLFCbyCompartment(listgr,namePCAcol="E1_compartment", padjVal=padjVal,
                                lfcVal=localLFC)

p1<-plotLFCbyCompartment(sigTbl,namePCAcol="E1_compartment", namePCA=paste(pcaSource,"E1"))

ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                "ABcomp_",pcaSource,"-E1_LFC_padj",
                                padjVal,"_lfc", localLFC,".pdf"),
                plot=p1, device="pdf",width=29,height=19,units="cm")

# test significance of LFC
for(s in unique(sigTbl$SMC)){
  print(s)
  sigTblss<-sigTbl[sigTbl$SMC==s,]
  resAOV<-aov(abs(log2FoldChange)~updown*E1_compartment,data=sigTblss)
  print(summary(resAOV))
  #  TukeyHSD(resAOV)
}


### second 366 eigenvector------

localLFC=0
sigTbl<-processLFCbyCompartment(listgr,namePCAcol="E2_compartment", padjVal=padjVal,
                                lfcVal=localLFC)

p1<-plotLFCbyCompartment(sigTbl,namePCAcol="E2_compartment", namePCA=paste(pcaSource,"E2"))

ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                "ABcomp_",pcaSource,"-E2_LFC_padj",
                                padjVal,"_lfc", localLFC,".pdf"),
                plot=p1, device="pdf",width=29,height=19,units="cm")



#############-
# Both eigen vectors combined------
#############-

localLFC=0
#localLFC=lfcVal

# unfiltered
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=localLFC,direction="both")
sigList<-lapply(sigList, "[", ,c("wormbaseID","E1_compartment","E2_compartment","log2FoldChange","padj"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList$E1.E2<-paste0(sigList$E1_compartment,"1.",sigList$E2_compartment,"2")
sigList<-sigList[!is.na(sigList$padj),]
sigList$E1.E2<-factor(sigList$E1.E2)
totalPerCat<-sigList %>% group_by(SMC,E1.E2) %>% summarise(count=n())

p1<-ggplot(totalPerCat,aes(x=E1.E2,y=count,group=SMC))+geom_bar(stat="identity")+
  facet_wrap(vars(SMC))

#yminmax=c(0,median(abs(sigList$log2FoldChange))+quantile(abs(sigList$log2FoldChange))[4]*2)
yminmax=c(0,1)
give.n <- function(x){
  return(c(y = 0, label = length(x)))
  # experiment with the multiplier to find the perfect position
}
p2<-ggplot(sigList,aes(x=E1.E2,y=abs(log2FoldChange))) +
  geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2),
               outlier.shape=NA, outlier.color="grey50")+
  facet_wrap(vars(SMC))+
  #ggtitle("Significantly changed genes by 366 E1.E2 compartments") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
  theme(legend.title = element_text(size=10)) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide="none") +
  stat_summary(fun.data = give.n, geom = "text",size=3,angle=90,
               position=position_dodge(width=0.6))




# upregulated
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=localLFC,direction="gt")

sigList<-lapply(sigList, "[", ,c("E1_compartment","E2_compartment","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList$E1.E2<-paste0(sigList$E1_compartment,"1.",sigList$E2_compartment,"2")

sigList<-sigList[!is.na(sigList$E1.E2),]
sigTbl<-sigList
sigTbl$updown<-"up"


# downregulated
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal, lfc= -localLFC, direction="lt")

sigList<-lapply(sigList, "[", ,c("E1_compartment","E2_compartment","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList$E1.E2<-paste0(sigList$E1_compartment,"1.",sigList$E2_compartment,"2")
sigList<-sigList[!is.na(sigList$E1.E2),]
sigList$updown<-"down"
sigTbl<-rbind(sigTbl,sigList)
sigTbl$E1.E2<-as.factor(sigTbl$E1.E2)
sigTbl$updown<-factor(sigTbl$updown, levels=c("up","down"))

yminmax=c(0,median(abs(sigTbl$log2FoldChange))+quantile(abs(sigTbl$log2FoldChange))[4]*2)
give.n <- function(x){
  return(c(y = 0, label = length(x)))
  # experiment with the multiplier to find the perfect position
}

p3<-ggplot(sigTbl,aes(x=E1.E2,y=abs(log2FoldChange),col=updown,fill=updown)) +
  geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2),outlier.shape=NA,
               outlier.color="grey50") +
  facet_wrap(facets=vars(SMC),nrow=2,strip.position="top") +
  ggtitle("Significantly changed genes by 366 E1.E2 compartments") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
  theme(legend.title = element_text(size=10)) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide="none") +
  stat_summary(fun.data = give.n, geom = "text",size=3,angle=90,
               position=position_dodge(width=0.6))

# test significance of LFC
summary(aov(abs(log2FoldChange)~updown+E1.E2,data=sigTbl))

p<-ggarrange(ggarrange(p1,p2,ncol=2,nrow=1),p3,nrow=1,ncol=1)
ggpubr::ggexport(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                "ABcomp_PMW366-E1.E2_LFC_padj",
                                padjVal,"_lfc", localLFC,".pdf"),
                plot=p, device="pdf",width=19,height=9,units="cm")





# AB compartments - sample specific ---------------------------------------

####
## sample specific compartments -----
####

RNAseqAndHiCsubset=c("aux_sdc3BG","dpy26","kle2","scc1","coh1")

if(all(RNAseqAndHiCsubset %in% useContrasts)){
  pcas<-data.frame(SMC=c("TEVonly","aux_sdc3BG","dpy26","kle2","scc1","coh1"),
                   strain =c("366","822","382","775","784","828"),
                   E1=NA, E2=NA)

  E1files=list.files(paste0(outPath,"/otherData"),
                            pattern="_merge_2000\\.oriented_E1\\.vecs\\.bw")
  E2files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.oriented_E2\\.vecs\\.bw")
  pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000.oriented_E1.vecs.bw")))]
  pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000.oriented_E2.vecs.bw")))]
  listgr<-NULL
  for (grp in RNAseqAndHiCsubset){
    #grp=useContrasts[1]
    salmon<-readRDS(file=paste0(outPath,"/rds/",fileNamePrefix,
                contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
    pca1<-import.bw(paste0(outPath,"/otherData/",pcas$E1[pcas$SMC==grp]))
    pca2<-import.bw(paste0(outPath,"/otherData/",pcas$E2[pcas$SMC==grp]))

    salmon<-salmon[!is.na(salmon$chr),]
    salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

    salmongr<-sort(salmongr)

    salmongr<-assignGRtoAB(salmongr,pca1,grName=grp,pcaName=paste0("E1"))
    salmongr<-assignGRtoAB(salmongr,pca2,grName=grp,pcaName=paste0("E2"))
    listgr[[grp]]<-salmongr
  }


  ####
  ## AB compartment by chromosome -----
  ####

  pcaSource="same"
  ### sample specific first eigenvector ------
  dfl<-processCountsPerChr(listgr,namePCAcol="E1_compartment")

  p1<-plotCountsPerChrPerCompartment(dfl, namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"))
  p1a<-plotFractionPerChrPerCompartment(dfl, namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"))

  dfl<-processUpDownCountsPerChrPerCompartment(listgr, namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"))

  p2<-plotUpDownByCompartment(dfl,namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"), compartment="A")
  p3<-plotUpDownByCompartment(dfl,namePCAcol="E1_compartment", namePCA=paste(pcaSource, "E1"), compartment="B")


  p<-ggpubr::ggarrange(p1,p1a,p2,p3,ncol=1,nrow=2)
  ggpubr::ggexport(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                                     "ABcomp_", pcaSource, "-E1_countsPerChr_padj",
                                     padjVal,"_lfc", lfcVal,".pdf"),
                   device="pdf",width=10,height=5, units="cm")


  ### sample specific second eigenvector ------
  dfl<-processCountsPerChr(listgr,namePCAcol="E2_compartment")

  p1<-plotCountsPerChrPerCompartment(dfl, namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"))
  p1a<-plotFractionPerChrPerCompartment(dfl, namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"))

  dfl<-processUpDownCountsPerChrPerCompartment(listgr, namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"))

  p2<-plotUpDownByCompartment(dfl,namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"), compartment="A")
  p3<-plotUpDownByCompartment(dfl,namePCAcol="E2_compartment", namePCA=paste(pcaSource, "E2"), compartment="B")

  p<-ggpubr::ggarrange(p1,p1a,p2,p3,ncol=1,nrow=2)
  ggpubr::ggexport(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                                     "ABcomp_",pcaSource,"-E2_countsPerChr_padj",
                                     padjVal,"_lfc", lfcVal,".pdf"),
                   device="pdf",width=10,height=5, units="cm")





  ####
  ## AB comp LFC -----
  ####

  localLFC=0
  pcaSource="same"
  ### sample specific fist eigenvector------
  sigTbl<-processLFCbyCompartment(listgr,namePCAcol="E1_compartment", padjVal=padjVal,
                                  lfcVal=localLFC)

  p1<-plotLFCbyCompartment(sigTbl,namePCAcol="E1_compartment", namePCA=paste(pcaSource,"E1"))

  ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                  "ABcomp_",pcaSource,"-E1_LFC_padj",
                                  padjVal,"_lfc", localLFC,".pdf"),
                  plot=p1, device="pdf",width=29,height=19,units="cm")

  # test significance of LFC
  for(s in unique(sigTbl$SMC)){
    print(s)
    sigTblss<-sigTbl[sigTbl$SMC==s,]
    resAOV<-aov(abs(log2FoldChange)~updown*E1_compartment,data=sigTblss)
    print(summary(resAOV))
    #  TukeyHSD(resAOV)
  }


  ###  sample specific second eigenvector------
  localLFC=0
  sigTbl<-processLFCbyCompartment(listgr,namePCAcol="E2_compartment", padjVal=padjVal,
                                  lfcVal=localLFC)

  p1<-plotLFCbyCompartment(sigTbl,namePCAcol="E2_compartment", namePCA=paste(pcaSource,"E2"))

  ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                  "ABcomp_",pcaSource,"-E2_LFC_padj",
                                  padjVal,"_lfc", localLFC,".pdf"),
                  plot=p1, device="pdf",width=29,height=19,units="cm")


  # AB compartments - switching ---------------------------------------------

  ####
  ## sample specific compartments - changes between TEVonly and cs
  ####
  pcas
  pcaSource="switch"
  listgr<-NULL
  for (grp in RNAseqAndHiCsubset){
    #grp=useContrasts[1]
    salmon<-readRDS(file=paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],
                                "_DESeq2_fullResults_p",padjVal,".rds"))


    pca1<-import.bw(paste0(outPath,"/otherData/",pcas$E1[pcas$SMC==grp]))
    pca2<-import.bw(paste0(outPath,"/otherData/",pcas$E2[pcas$SMC==grp]))
    pca1wt<-import.bw(paste0(outPath,"/otherData/",pcas$E1[pcas$SMC=="TEVonly"]))
    pca2wt<-import.bw(paste0(outPath,"/otherData/",pcas$E2[pcas$SMC=="TEVonly"]))

    salmon<-salmon[!is.na(salmon$chr),]
    salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

    salmongr<-sort(salmongr)

    salmongr<-assignGRtoAB(salmongr,pca1,grName=grp,pcaName=paste0("E1"))
    salmongr<-assignGRtoAB(salmongr,pca2,grName=grp,pcaName=paste0("E2"))
    salmongr<-assignGRtoAB(salmongr,pca1wt,grName=grp,pcaName=paste0("E1wt"))
    salmongr<-assignGRtoAB(salmongr,pca2wt,grName=grp,pcaName=paste0("E2wt"))

    salmongr$E1_switch<-factor(paste0(salmongr$E1wt_compartment,salmongr$E1_compartment),levels=c("AA","BB","AB","BA"))

    salmongr$E2_switch<-factor(paste0(salmongr$E2wt_compartment,salmongr$E2_compartment),levels=c("AA","BB","AB","BA"))

    listgr[[grp]]<-salmongr
  }



  ################-
  ## AB compartment by chromosome - switching between TEVonly and cs -----
  #################-


  pcaSource="switch"

  ### switch first Eigenvector ------
  dfl<-processCountsPerChr(listgr,namePCAcol="E1_switch")

  p1<-plotCountsPerChrPerCompartment(dfl, namePCAcol="E1_switch", namePCA=paste(pcaSource, "E1"))
  p1a<-plotFractionPerChrPerCompartment(dfl, namePCAcol="E1_switch", namePCA=paste(pcaSource, "E1"))

  dfl<-processUpDownCountsPerChrPerCompartment(listgr, namePCAcol="E1_switch", namePCA=paste(pcaSource, "E1"))

  p2<-plotUpDownByCompartment(dfl,namePCAcol="E1_switch", namePCA=paste(pcaSource, "E1"), compartment="AA")
  p3<-plotUpDownByCompartment(dfl,namePCAcol="E1_switch", namePCA=paste(pcaSource, "E1"), compartment="BB")
  p4<-plotUpDownByCompartment(dfl,namePCAcol="E1_switch", namePCA=paste(pcaSource, "E1"), compartment="AB")
  p5<-plotUpDownByCompartment(dfl,namePCAcol="E1_switch", namePCA=paste(pcaSource, "E1"), compartment="BA")

  p<-ggpubr::ggarrange(p1,p1a,p2,p3,p4,p5,ncol=1,nrow=2)
  ggpubr::ggexport(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                                     "ABcomp_", pcaSource, "-E1_countsPerChr_padj",
                                     padjVal,"_lfc", lfcVal,".pdf"),
                   device="pdf",width=10,height=5, units="cm")


  ### switch second eigenvector ------
  dfl<-processCountsPerChr(listgr,namePCAcol="E2_switch")

  p1<-plotCountsPerChrPerCompartment(dfl, namePCAcol="E2_switch", namePCA=paste(pcaSource, "E2"))
  p1a<-plotFractionPerChrPerCompartment(dfl, namePCAcol="E2_switch", namePCA=paste(pcaSource, "E2"))

  dfl<-processUpDownCountsPerChrPerCompartment(listgr, namePCAcol="E2_switch", namePCA=paste(pcaSource, "E2"))

  p2<-plotUpDownByCompartment(dfl,namePCAcol="E2_switch", namePCA=paste(pcaSource, "E2"), compartment="AA")
  p3<-plotUpDownByCompartment(dfl,namePCAcol="E2_switch", namePCA=paste(pcaSource, "E2"), compartment="BB")
  p4<-plotUpDownByCompartment(dfl,namePCAcol="E2_switch", namePCA=paste(pcaSource, "E2"), compartment="AB")
  p5<-plotUpDownByCompartment(dfl,namePCAcol="E2_switch", namePCA=paste(pcaSource, "E2"), compartment="BA")

  p<-ggpubr::ggarrange(p1,p1a,p2,p3,p4,p5,ncol=1,nrow=2)
  ggpubr::ggexport(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                                     "ABcomp_",pcaSource,"-E2_countsPerChr_padj",
                                     padjVal,"_lfc", lfcVal,".pdf"),
                   device="pdf",width=10,height=5, units="cm")



  ####
  ## AB comp LFC - switching between TEVonly and cs -----
  ####

  localLFC=0
  pcaSource="switch"
  ### sample specific fist eigenvector------
  sigTbl<-processLFCbyCompartment(listgr,namePCAcol="E1_switch", padjVal=padjVal,
                                  lfcVal=localLFC)

  p1<-plotLFCbyCompartment(sigTbl,namePCAcol="E1_switch", namePCA=paste(pcaSource,"E1"))

  ggpubr::ggexport(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                  "ABcomp_",pcaSource,"-E1_LFC_padj",
                                  padjVal,"_lfc", localLFC,".pdf"),
                  plot=p1, device="pdf",width=10,height=5,units="cm")

  # test significance of LFC
  for(s in unique(sigTbl$SMC)){
    print(s)
    sigTblss<-sigTbl[sigTbl$SMC==s,]
    resAOV<-aov(abs(log2FoldChange)~updown*E1_switch,data=sigTblss)
    print(summary(resAOV))
    #  TukeyHSD(resAOV)
  }


  ###  sample specific second eigenvector------
  localLFC=0
  sigTbl<-processLFCbyCompartment(listgr,namePCAcol="E2_switch", padjVal=padjVal,
                                  lfcVal=localLFC)

  p1<-plotLFCbyCompartment(sigTbl,namePCAcol="E2_switch", namePCA=paste(pcaSource,"E2"))

  ggpubr::ggexport(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                  "ABcomp_",pcaSource,"-E2_LFC_padj",
                                  padjVal,"_lfc", localLFC,".pdf"),
                  plot=p1, device="pdf",width=10,height=5,units="cm")
}



##################-
## 366 TPM in AB compartments of different HiCs -----
#################-
RNAseqAndHiCsubset=c("aux_sdc3BG","dpy26","kle2","scc1","coh1")

if(all(RNAseqAndHiCsubset %in% useContrasts)){
  pcas<-data.frame(SMC=c("TEVonly","aux_sdc3BG","dpy26","kle2","scc1","coh1"),
                   strain =c("366","822","382","775","784","828"),
                   E1=NA, E2=NA)

  tpm366<-import(paste0(outPath,"/tracks/PMW366_TPM_avr.bedgraph"),format="bedgraph")
  cov366<-coverage(tpm366,weight="score")

  E1files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.oriented_E1\\.vecs\\.bw")
  E2files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.oriented_E2\\.vecs\\.bw")
  pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.oriented_E1\\.vecs\\.bw")))]
  pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.oriented_E2\\.vecs\\.bw")))]
  listdf<-NULL
  for (grp in pcas$SMC){
    #grp=useContrasts[1]
    # salmon<-import.bed(file=paste0(outPath,"/rds/",fileNamePrefix,
    #                             contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
    pca1<-import.bw(paste0(outPath,"/otherData/",pcas$E1[pcas$SMC==grp]))
    pca2<-import.bw(paste0(outPath,"/otherData/",pcas$E2[pcas$SMC==grp]))

    seqlevels(pca1)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
    seqlevels(pca2)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)

    pca1<-binnedAverage(pca1,cov366,varname="tpm366")
    pca2<-binnedAverage(pca2,cov366,varname="tpm366")

    df1<-data.frame(pca1)
    df2<-data.frame(pca2)

    df1$pca<-"E1"
    df2$pca<-"E2"

    df1$compartment<-NA
    df2$compartment<-NA
    df1$compartment<-ifelse(df1$score>0,"A","B")
    df2$compartment<-ifelse(df2$score>0,"A","B")

    df<-rbind(df1,df2)
    df$compartment<-factor(df$compartment)
    df$SMC<-grp

    listdf[[grp]]<-df
  }

  df<-do.call(rbind,listdf)
  df$SMC<-factor(df$SMC,levels=pcas$SMC)

  p<-ggplot(df,aes(x=compartment,y=log2(tpm366))) +
    geom_boxplot(outlier.shape=NA) + facet_grid(cols=vars(SMC),rows=vars(pca)) +
    ylim(c(-15,15)) + theme_bw()+ geom_hline(yintercept=0,col="red")+
    ggtitle(paste0("366 TPM in different bins of pca"))

  ggsave(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                           "eigenValAll_366tpm.pdf"),
         device="pdf",width=29,height=19, units="cm")

  corMethod="pearson"
  tpmThresh=0
  allBins<-nrow(df)
  df<-df[df$tpm366>=tpmThresh,]
  fracKept<-nrow(df)/allBins # 0.79 of bins
  p<-ggplot(df,aes(x=score,y=log2(tpm366))) +
    #geom_point(colour="#00007733") +
    #geom_density2d_filled()+
    #geom_hex(bins=100)+
    geom_bin2d(bins=100)+
    #stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    facet_grid(rows=vars(pca),cols=vars(SMC)) + geom_smooth(method="lm") +
    stat_cor(label.x = -1.5, label.y = 18, size=3,method=corMethod) +
    ggtitle(paste0(corMethod," correlation of PCA eigen value vs PMW366 TPM (for bins > ",
                   formatC(tpmThresh,big.mark=",",format="G"),"tpm, ",
                   round(100*fracKept,1),"% of bins)")) +
    ylim(c(-20,20)) +xlim(c(-1.5,1.5))
  ggsave(p,filename=paste0(outPath, "/plots/",outputNamePrefix,corMethod,
                          "Cor_eigenValAll_366tpm",tpmThresh,".pdf"),
         device="pdf",width=29,height=14, units="cm")
}



##################-
## sampleSpecific TPM in AB compartments of different HiCs -----
#################-
RNAseqAndHiCsubset=c("aux_sdc3BG","dpy26","kle2","scc1","coh1")

if(all(RNAseqAndHiCsubset %in% useContrasts)){
  pcas<-data.frame(SMC=c("TEVonly","aux_sdc3BG","dpy26","kle2","scc1","coh1"),
                   strain =c("366","822","382","775","784","828"),
                   E1=NA, E2=NA)


  E1files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.oriented_E1\\.vecs\\.bw")
  E2files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.oriented_E2\\.vecs\\.bw")
  pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.oriented_E1\\.vecs\\.bw")))]
  pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.oriented_E2\\.vecs\\.bw")))]
  listdf<-NULL
  g=1
  for (g in 1:nrow(pcas)){
    grp<-pcas$SMC[g]
    #grp=useContrasts[1]
    # salmon<-import.bed(file=paste0(outPath,"/rds/",fileNamePrefix,
    #                             contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

    tpmFile<-paste0(outPath,"/tracks/PMW",pcas$strain[g],"_TPM_avr.bedgraph")
    print(tpmFile)
    tpm<-import(tpmFile,format="bedgraph")
    covstrain<-coverage(tpm,weight="score")

    pca1<-import.bw(paste0(outPath,"/otherData/",pcas$E1[pcas$SMC==grp]))
    pca2<-import.bw(paste0(outPath,"/otherData/",pcas$E2[pcas$SMC==grp]))

    seqlevels(pca1)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
    seqlevels(pca2)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)

    pca1<-binnedAverage(pca1,covstrain,varname="tpm")
    pca2<-binnedAverage(pca2,covstrain,varname="tpm")

    df1<-data.frame(pca1)
    df2<-data.frame(pca2)

    df1$pca<-"E1"
    df2$pca<-"E2"

    df1$compartment<-NA
    df2$compartment<-NA
    df1$compartment<-ifelse(df1$score>0,"A","B")
    df2$compartment<-ifelse(df2$score>0,"A","B")

    df<-rbind(df1,df2)
    df$compartment<-factor(df$compartment)
    df$SMC<-grp

    listdf[[grp]]<-df
  }

  df<-do.call(rbind,listdf)
  df$SMC<-factor(df$SMC,levels=pcas$SMC)

  p<-ggplot(df,aes(x=compartment,y=log2(tpm))) +
    geom_boxplot(outlier.shape=NA) + facet_grid(SMC~pca) +
    ylim(c(-15,15)) + theme_bw()+ geom_hline(yintercept=0,col="red")+
    ggtitle(paste0("TPM in different bins of pca"))

  ggsave(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                           "eigenValAll_sameTPM.pdf"),
         device="pdf",width=29,height=19, units="cm")

  corMethod="pearson"
  tpmThresh=0
  allBinNum<-df %>% dplyr::group_by(SMC) %>% dplyr::summarise(count=n())
  allBinNum$autosomal<-df %>% dplyr::group_by(SMC) %>% filter(seqnames!="chrX") %>% dplyr::summarise(count=n())
  df<-df[df$tpm>=tpmThresh,]
  keptBinNum<-df %>% dplyr::group_by(SMC) %>% dplyr::summarise(count=n())
  keptBinNum$percent<-round(100*keptBinNum$count/allBinNum$count,1)
  p<-ggplot(df,aes(x=score,y=log2(tpm))) +
    #geom_point(colour="#00007733") +
    #geom_density2d_filled()+
    #geom_hex(bins=100)+
    geom_bin2d(bins=100)+
    #stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    facet_grid(rows=vars(pca),cols=vars(SMC)) + geom_smooth(method="lm") +
    stat_cor(label.x = -1.5, label.y = 18, size=3,method=corMethod) +
    ggtitle(paste0(corMethod," correlation of PCA eigen value vs TPM (for bins > ",
                   formatC(tpmThresh,big.mark=",",format="G"),"tpm, >",
                   min(keptBinNum$percent),"% of bins)")) +
    ylim(c(-20,20)) +xlim(c(-1.5,1.5))
  ggsave(p,filename=paste0(outPath, "/plots/",outputNamePrefix,corMethod,
                           "Cor_eigenValAll_sameTPM",tpmThresh,".pdf"),
         device="pdf",width=29,height=14, units="cm")

  keptBinNum<-df[df$seqnames!="chrX",]  %>% dplyr::group_by(SMC) %>% dplyr::summarise(countChrA=n()) %>%mutate(percentChrA=round(100*countChrA/allBinNum$autosomal$count,1))
  p1<-ggplot(df[df$seqnames!="chrX",],aes(x=score,y=log2(tpm))) +
    #geom_point(colour="#00007733") +
    #geom_density2d_filled()+
    #geom_hex(bins=100)+
    geom_bin2d(bins=100)+
    #stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    facet_grid(rows=vars(pca),cols=vars(SMC)) + geom_smooth(method="lm") +
    stat_cor(label.x = -1.5, label.y = 18, size=3,method=corMethod) +
    ggtitle(paste0(corMethod," correlation of autosomal PCA eigen value vs TPM (for bins > ",
                   formatC(tpmThresh,big.mark=",",format="G"),"tpm, >",
                   min(keptBinNum$percentChrA),"% of chrA bins)")) +
    ylim(c(-20,20)) +xlim(c(-1.5,1.5))

  ggsave(p1,filename=paste0(outPath, "/plots/",outputNamePrefix,corMethod,
                           "Cor_eigenValAll_sameTPM",tpmThresh,"_chrA.pdf"),
         device="pdf",width=29,height=14, units="cm")

}





###########################-
# compartments - digitized ----------------------------------------------------
###########################-

## 366TPM in digitized compartments of different HiCs -----

### Autosomes ------
RNAseqAndHiCsubset=c("aux_sdc3BG","dpy26","kle2","scc1","coh1")

if(all(RNAseqAndHiCsubset %in% useContrasts)){
  pcas<-data.frame(SMC=c("wt","TEVonly","aux_sdc3BG","dpy26","kle2","scc1","coh1"),
                   strain =c("N2","366","822","382","775","784","828"),
                   E1=NA, E2=NA)

  tpm366<-import(paste0(outPath,"/tracks/PMW366_TPM_avr.bedgraph"),format="bedgraph")
  cov366<-coverage(tpm366,weight="score")

  E1files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.saddle_trans_A_noX_E1\\.digitized\\.rds")
  E2files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.saddle_trans_A_noX_E2\\.digitized\\.rds")
  pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.saddle_trans_A_noX_E1\\.digitized\\.rds")))]
  pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.saddle_trans_A_noX_E2\\.digitized\\.rds")))]
  listdf<-NULL
  for (grp in pcas$SMC){
    pca1<-readRDS(paste0(outPath,"/otherData/",pcas$E1[pcas$SMC==grp]))
    pca2<-readRDS(paste0(outPath,"/otherData/",pcas$E2[pcas$SMC==grp]))

    pca1<-binnedAverage(pca1,cov366,varname="tpm366")
    pca2<-binnedAverage(pca2,cov366,varname="tpm366")

    df1<-data.frame(pca1)
    df2<-data.frame(pca2)

    df<-rbind(df1,df2)
    #df$compartment<-factor(df$compartment)
    df$SMC<-grp

    listdf[[grp]]<-df
  }

  df<-do.call(rbind,listdf)
  df$SMC<-factor(df$SMC,levels=pcas$SMC)
  df$bin<-factor(df$bin,levels=1:50)


  p<-ggplot(df,aes(x=bin,y=log2(tpm366))) +
    geom_boxplot(outlier.shape=NA) + facet_grid(SMC~pca) +
    ylim(c(-15,15)) + theme_bw()+ geom_hline(yintercept=0,col="red")+
    ggtitle(paste0("366 TPM in different ausotomal bins of digitized pca"))

  ggsave(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                  "digitizedCompAll_chrA_366tpm.pdf"),
  device="pdf",width=29,height=19, units="cm")

  subdf<-df[df$SMC %in% c("wt","TEVonly"),]
  subdf<-subdf[subdf$bin %in% 1:50,]
  p<-ggplot(subdf,aes(x=bin,y=log2(tpm366),fill=bin)) +
    geom_boxplot(outlier.shape=NA,size=0.1,fill="lightblue") + facet_grid(SMC~pca)+
    ylim(c(-12,12)) + theme_bw()+
    #scale_fill_manual(values=c("white","grey70"))+
    geom_hline(yintercept=0,col="red")+
    ggtitle(paste0("366 TPM in different autosomal bins of digitized pca")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank())
  ggsave(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                           "digitizedCompControls_chrA_366tpm.pdf"),
         device="pdf",width=10,height=12, units="cm")

}


### chrX -------
RNAseqAndHiCsubset=c("aux_sdc3BG","dpy26","kle2","scc1","coh1")

if(all(RNAseqAndHiCsubset %in% useContrasts)){
  pcas<-data.frame(SMC=c("wt","TEVonly","aux_sdc3BG","dpy26","kle2","scc1","coh1"),
                   strain =c("N2","366","822","382","775","784","828"),
                   E1=NA, E2=NA)

  tpm366<-import(paste0(outPath,"/tracks/PMW366_TPM_avr.bedgraph"),format="bedgraph")
  cov366<-coverage(tpm366,weight="score")

  E1files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.saddle_cis_X_noA_E1\\.digitized\\.rds")
  E2files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.saddle_cis_X_noA_E2\\.digitized\\.rds")
  pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.saddle_cis_X_noA_E1\\.digitized\\.rds")))]
  pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.saddle_cis_X_noA_E2\\.digitized\\.rds")))]
  listdf<-NULL
  for (grp in pcas$SMC){
    pca1<-readRDS(paste0(outPath,"/otherData/",pcas$E1[pcas$SMC==grp]))
    pca2<-readRDS(paste0(outPath,"/otherData/",pcas$E2[pcas$SMC==grp]))

    pca1<-binnedAverage(pca1,cov366,varname="tpm366")
    pca2<-binnedAverage(pca2,cov366,varname="tpm366")

    df1<-data.frame(pca1)
    df2<-data.frame(pca2)

    df<-rbind(df1,df2)
    #df$compartment<-factor(df$compartment)
    df$SMC<-grp

    listdf[[grp]]<-df
  }

  df<-do.call(rbind,listdf)
  df$SMC<-factor(df$SMC,levels=pcas$SMC)
  df$bin<-factor(df$bin,levels=1:50)

  p<-ggplot(df,aes(x=bin,y=log2(tpm366))) +
    geom_boxplot(outlier.shape=NA) + facet_grid(SMC~pca) +
    ylim(c(-15,15)) + theme_bw()+ geom_hline(yintercept=0,col="red")+
    ggtitle(paste0("366 TPM in different chrX bins of digitized pca"))

  ggsave(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                           "digitizedCompAll_chrX_366tpm.pdf"),
         device="pdf",width=29,height=19, units="cm")

  subdf<-df[df$SMC %in% c("wt","TEVonly"),]
  subdf<-subdf[subdf$bin %in% 1:50,]
  p<-ggplot(subdf,aes(x=bin,y=log2(tpm366),fill=bin)) +
    geom_boxplot(outlier.shape=NA,size=0.1,fill="lightblue") + facet_grid(SMC~pca)+
    ylim(c(-12,12)) + theme_bw()+
    #scale_fill_manual(values=c("white","grey70"))+
    geom_hline(yintercept=0,col="red")+
    ggtitle(paste0("366 TPM in different chrX bins of digitized pca")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank())
  ggsave(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                           "digitizedCompControls_chrX_366tpm.pdf"),
         device="pdf",width=10,height=12, units="cm")

}


## log baseMean in digitized compartments of different HiCs ------

RNAseqAndHiCsubset=c("aux_sdc3BG","dpy26","kle2","scc1","coh1")

if(all(RNAseqAndHiCsubset %in% useContrasts)){
  pcas<-data.frame(SMC=c("TEVonly","aux_sdc3BG","dpy26","kle2","scc1","coh1"),
                   strain =c("366","822","382","775","784","828"),
                   E1=NA, E2=NA)

  E1files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.saddle_trans_A_noX_E1\\.digitized\\.rds")
  E2files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.saddle_trans_A_noX_E2\\.digitized\\.rds")
  pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.saddle_trans_A_noX_E1\\.digitized\\.rds")))]
  pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.saddle_trans_A_noX_E2\\.digitized\\.rds")))]


  salmon<-readRDS(file=paste0(paste0(outPath,"/rds/",fileNamePrefix,
                                     contrastNames[["dpy26"]],"_DESeq2_fullResults_p",padjVal,".rds")))

  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

  baseMean<-coverage(salmongr,weight="baseMean")

  listdf<-NULL
  for (grp in RNAseqAndHiCsubset){
    #grp=useContrasts[1]
    pca1<-readRDS(paste0(outPath,"/otherData/",pcas$E1[pcas$SMC==grp]))
    pca2<-readRDS(paste0(outPath,"/otherData/",pcas$E2[pcas$SMC==grp]))

    pca1<-binnedAverage(pca1,cov366,varname="baseMean")
    pca2<-binnedAverage(pca2,cov366,varname="baseMean")

    df1<-data.frame(pca1)
    df2<-data.frame(pca2)

    df<-rbind(df1,df2)
    df$bin<-factor(df$bin,levels=1:50)
    df$SMC<-grp

    listdf[[grp]]<-df
  }

  df<-do.call(rbind,listdf)
  df$SMC<-factor(df$SMC,levels=pcas$SMC)

  p<-ggplot(df,aes(x=bin,y=baseMean)) +
    geom_boxplot(outlier.shape=NA) + facet_grid(SMC~pca) +
    ylim(c(0,60)) + theme_bw()+ #geom_hline(yintercept=0,col="red")+
    ggtitle(paste0("baseMean RNAseq in bins of different strain digitized pca"))

  ggsave(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                           "digitizedCompAll_baseMean.pdf"),
         device="pdf",width=29,height=19, units="cm")
}

######### lfc in digitized compartments of different HiCs------

RNAseqAndHiCsubset=c("aux_sdc3BG","dpy26","kle2","scc1","coh1")

if(all(RNAseqAndHiCsubset %in% useContrasts)){
  pcas<-data.frame(SMC=c("TEVonly","aux_sdc3BG","dpy26","kle2","scc1","coh1"),
                   strain =c("366","822","382","775","784","828"),
                   E1=NA, E2=NA)

  E1files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.saddle_trans_E1\\.digitized\\.rds")
  E2files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.saddle_trans_E2\\.digitized\\.rds")
  pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.saddle_trans_E1\\.digitized\\.rds")))]
  pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.saddle_trans_E2\\.digitized\\.rds")))]
  listdf<-NULL
  for (grp in RNAseqAndHiCsubset){
    #grp=useContrasts[1]
    salmon<-readRDS(file=paste0(paste0(outPath,"/rds/",fileNamePrefix,
                                       contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))

    salmon<-salmon[!is.na(salmon$chr),]
    salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

    lfc<-coverage(salmongr,weight="log2FoldChange")

    pca1<-readRDS(paste0(outPath,"/otherData/",pcas$E1[pcas$SMC==grp]))
    pca2<-readRDS(paste0(outPath,"/otherData/",pcas$E2[pcas$SMC==grp]))

    pca1<-binnedAverage(pca1,cov366,varname="lfc")
    pca2<-binnedAverage(pca2,cov366,varname="lfc")

    df1<-data.frame(pca1)
    df2<-data.frame(pca2)

    df<-rbind(df1,df2)
    df$pca<-factor(df$pca)
    df$SMC<-grp

    listdf[[grp]]<-df
  }

  df<-do.call(rbind,listdf)
  df$SMC<-factor(df$SMC,levels=pcas$SMC)
  df$bin<-factor(df$bin,levels=1:50)
  p<-ggplot(df,aes(x=bin,y=lfc)) +
    geom_boxplot(outlier.shape=NA) + facet_grid(SMC~pca) +
    ylim(c(0,60)) + theme_bw()+ #geom_hline(yintercept=0,col="red")+
    ggtitle(paste0("Log2FC in different bins of digitized pca (lfc/pca same strain)"))
  ggsave(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                           "digitizedCompAll_lfcAll.pdf"),
         device="pdf",width=29,height=19, units="cm")
}



######### lfc in digitized compartments of 366 HiC-----

RNAseqAndHiCsubset=c("aux_sdc3BG","dpy26","kle2","scc1","coh1")

if(all(RNAseqAndHiCsubset %in% useContrasts)){
  pcas<-data.frame(SMC=c("TEVonly","aux_sdc3BG","dpy26","kle2","scc1","coh1"),
                   strain =c("366","822","382","775","784","828"),
                   E1=NA, E2=NA)


  E1files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.saddle_trans_E1\\.digitized\\.rds")
  E2files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.saddle_trans_E2\\.digitized\\.rds")
  pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.saddle_trans_E1\\.digitized\\.rds")))]
  pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.saddle_trans_E2\\.digitized\\.rds")))]
  listdf<-NULL
  for (grp in RNAseqAndHiCsubset){
    #grp=useContrasts[1]
    salmon<-readRDS(file=paste0(paste0(outPath,"/rds/",fileNamePrefix,
                                       contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))

    salmon<-salmon[!is.na(salmon$chr),]
    salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

    lfc<-coverage(salmongr,weight="log2FoldChange")

    pca1<-readRDS(paste0(outPath,"/otherData/",pcas$E1[pcas$SMC=="TEVonly"]))
    pca2<-readRDS(paste0(outPath,"/otherData/",pcas$E2[pcas$SMC=="TEVonly"]))

    pca1<-binnedAverage(pca1,cov366,varname="lfc")
    pca2<-binnedAverage(pca2,cov366,varname="lfc")

    df1<-data.frame(pca1)
    df2<-data.frame(pca2)

    df<-rbind(df1,df2)
    df$pca<-factor(df$pca)
    df$SMC<-grp

    listdf[[grp]]<-df
  }

  df<-do.call(rbind,listdf)
  df$SMC<-factor(df$SMC,levels=pcas$SMC)

  p<-ggplot(df,aes(x=bin,y=lfc)) +
    geom_boxplot(outlier.shape=NA) + facet_grid(SMC~pca) +
    ylim(c(0,60)) + theme_bw()+ #geom_hline(yintercept=0,col="red")+
    ggtitle(paste0("Log2FC of different strain in bins of TEVonly digitized pca"))
  ggsave(p,filename=paste0(outPath, "/plots/",outputNamePrefix,
                           "digitizedComp366_lfcAll.pdf"),
         device="pdf",width=29,height=19, units="cm")
}



# seqplots heatmaps and averages ---------------------------------------------


#####################


# ######################-
# # TADs ----
# ######################-
# # loops<-rtracklayer::import(paste0(outPath,"/otherData/N2.allValidPairs.hic.5-10kbLoops.bedpe"), format="bedpe")
# # head(loops)
# # #extract the separate anchors
# # grl<-zipup(loops)
# # anchor1<-unlist(grl)[seq(1,2*length(grl),2)]
# # anchor2<-unlist(grl)[seq(2,2*length(grl),2)]
#
# # new SIP loops
# loops<-read.delim(paste0(outPath,"/otherData/10kbLoops.txt"),header=T)
# anchor1<-GRanges(seqnames=loops$chromosome1,ranges=IRanges(start=loops$x1,end=loops$x2-1),
#                  strand="+")
# anchor2<-GRanges(seqnames=loops$chromosome2,ranges=IRanges(start=loops$y1,end=loops$y2-1),
#                  strand="-")
# mcols(anchor1)<-loops[,c(8:15)]
# mcols(anchor2)<-loops[,c(8:15)]
# anchor1$loopNum<-paste0("loop",1:length(anchor1))
# anchor2$loopNum<-paste0("loop",1:length(anchor2))
#
# #make TADs
# tads_in<-GRanges(seqnames=seqnames(anchor1),IRanges(start=end(anchor1)+1,end=start(anchor2)-1))
#
# tads_in<-tads_in[width(tads_in)>100000]
#
# #tads<-reduce(tads)
# xtads<-tads_in[seqnames(tads_in)=="chrX"]
# atads<-tads_in[seqnames(tads_in)!="chrX"]
#
# flankSize<-200000
#
# smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
#                   useContrasts,"_lfc.bw")
# names(smcRNAseq)<-useContrasts
#
# pdf(file=paste0(outPath,"/plots/",outputNamePrefix,"tads-chrX_flank",
#                   flankSize/1000,"kb.pdf"), width=11,
#       height=9, paper="a4r")
#
#
# p<-getPlotSetArray(tracks=c(smcRNAseq),
#                    features=c(xtads),
#                    refgenome="ce11", bin=10000L, xmin=flankSize,
#                    xmax=flankSize, type="af",
#                    xanchored=median(width(xtads)))
#
#
#
# dd<-plotHeatmap(p,plotz=F)
# heatmapQuantiles<-sapply(dd$HLST,quantile,c(0.05,0.95),na.rm=T)
# roworder<-rev(order(lapply(dd$HLST,rowSums,na.rm=T)$X1))
# minVal<-min(heatmapQuantiles[1,])
# maxVal<-max(heatmapQuantiles[2,])
# #layout(matrix(c(1), nrow = 2, ncol = 1, byrow = TRUE))
# plotAverage(p,ylim=c(-0.5,1),main="ChrX TADs",
#             error.estimates=ifelse(length(useContrasts>3),F,T))
# plotHeatmap(p,main="ChrX TADs", plotScale="no", sortrows=T,
#             clusters=1L,autoscale=F,zmin=minVal, zmax=maxVal,
#             indi=F, sort_mids=T,sort_by=c(T,F,F),
#             clspace=c("#00008B", "#FFFFE0","#8B0000"))
#
#
# # reduced tads
# p<-getPlotSetArray(tracks=c(smcRNAseq),
#                    features=c(reduce(xtads)),
#                    refgenome="ce11", bin=10000L, xmin=flankSize,
#                    xmax=flankSize, type="af",
#                    xanchored=median(width(xtads)))
#
#
# dd<-plotHeatmap(p,plotz=F)
# heatmapQuantiles<-sapply(dd$HLST,quantile,c(0.05,0.95),na.rm=T)
# roworder<-rev(order(lapply(dd$HLST,rowSums,na.rm=T)$X1))
# minVal<-min(heatmapQuantiles[1,])
# maxVal<-max(heatmapQuantiles[2,])
# #layout(matrix(c(1), nrow = 2, ncol = 1, byrow = TRUE))
# plotAverage(p,ylim=c(-0.5,1),main="ChrX TADs (reduced)",error.estimates=F)
# plotHeatmap(p,main="ChrX TADs (reduced)", plotScale="no", sortrows=T,
#             clusters=1L,autoscale=F,zmin=minVal, zmax=maxVal,
#             indi=F, sort_mids=T,sort_by=c(T,F,F),
#             clspace=c("#00008B", "#FFFFE0","#8B0000"))
#
#
# dev.off()




####
## loops-----
####

# smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
#                   useContrasts,"_lfc.bw")
#
# loops<-rtracklayer::import(paste0(outPath,"/otherData/N2.allValidPairs.hic.5-10kbLoops.bedpe"),
#                            format="bedpe")
# head(loops)
# #extract the separate anchors
# grl<-zipup(loops)
# anchor1<-unlist(grl)[seq(1,2*length(grl),2)]
# anchor2<-unlist(grl)[seq(2,2*length(grl),2)]
#
#
#
# anchors<-reduce(sort(c(anchor1,anchor2)))
# olap_gr<-anchors
# target_size<-max(width(anchors))
# window_size<-target_size/2
# target_size=round(target_size/window_size)*window_size
# olap_gr<-resize(olap_gr,width=target_size,fix="center")
# bw_gr <- ssvFetchBigwig(smcRNAseq[3], olap_gr, win_size = 100,
#                        unique_names=useContrasts[3],win_method="summary")
# bw_gr$chr<-factor(gsub("chr","",as.vector(seqnames(bw_gr))))
# ssvSignalHeatmap(bw_gr,fill_limits=c(-1,1),
#                 perform_clustering="no",cluster_="chr",
#                 within_order_strategy="sort")
#
# library(ComplexHeatmap)
# library(seqsetvis)
#
# grpbw<-import.bw(smcRNAseq[3])
# bwsig<-viewGRangesWinSummary_dt(grpbw,olap_gr,n_tiles=100)
# mat<-matrix(bwsig$y,ncol=100,byrow=T)
# col_fun = circlize::colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
#
# p1<-Heatmap(mat,cluster_columns=F,cluster_rows=F,
#             row_split=bwsig$seqnames[seq(1,dim(bwsig)[1],100)],
#             show_row_dend = F)#, col=col_fun(c(seq(-1,1,1))),na_col="grey")
# p1
#
# p1+p2


# # ## old loops
# # loops<-rtracklayer::import(paste0(outPath,"/otherData/N2.allValidPairs.hic.5-10kbLoops.bedpe"),
# #                            format="bedpe")
# # head(loops)
# # #extract the separate anchors
# # grl<-zipup(loops)
# # anchor1<-unlist(grl)[seq(1,2*length(grl),2)]
# # anchor2<-unlist(grl)[seq(2,2*length(grl),2)]
#
#
# ### new SIP loops
# loops<-read.delim(paste0(outPath,"/otherData/10kbLoops.txt"),header=T)
# anchor1<-GRanges(seqnames=loops$chromosome1,ranges=IRanges(start=loops$x1,end=loops$x2-1))
# anchor2<-GRanges(seqnames=loops$chromosome2,ranges=IRanges(start=loops$y1,end=loops$y2-1))
# mcols(anchor1)<-loops[,c(8:15)]
# mcols(anchor2)<-loops[,c(8:15)]
# anchor1$loopNum<-paste0("loop",1:length(anchor1))
# anchor2$loopNum<-paste0("loop",1:length(anchor2))
#
# # anchors<-c(anchor1,anchor2)
# # #anchors<-reduce(anchors,min.gapwidth=0L)
# # seqlevels(anchors)<-seqlevels(Celegans)[1:6]
#
#
# #make TADs
# tads<-GRanges(seqnames=seqnames(anchor1),IRanges(start=start(anchor1),end=end(anchor2)))
# head(tads)
# sort(width(tads))
# tads<-reduce(tads)
# sort(width(tads))
#
# # find regions not in tads
# notads<-gaps(tads)
# sort(width(notads))
#
#
# plotList<-list()
# #grp=useContrasts[3]
# for (grp in useContrasts){
#   salmon<-readRDS(file=paste0(paste0(outPath,"/rds/",fileNamePrefix,
#                                      contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
#
#   salmon<-salmon[!is.na(salmon$chr),]
#   salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)
#
#   salmongr<-sort(salmongr)
#
#   ol<-findOverlaps(salmongr,tads,type="within")
#   genesInTads<-salmongr[queryHits(ol)]
#
#   ol<-findOverlaps(salmongr,notads)
#   genesNotTads<-salmongr[queryHits(ol)]
#
#   genesInTads$TADs<-"inside"
#   genesNotTads$TADs<-"outside"
#   df<-data.frame(c(genesInTads,genesNotTads))
#   df<-df%>%dplyr::group_by(seqnames,TADs)%>%dplyr::mutate(count=n())
#
#   plotList[[grp]]<-ggplot(df,aes(x=TADs,y=log2FoldChange,fill=TADs))+
#     geom_boxplot(notch=T,outlier.shape=NA,varwidth=T)+
#     facet_grid(.~seqnames) +ylim(c(-1,1))+
#     ggtitle(grp)
# }
# p<-gridExtra::marrangeGrob(plotList,ncol=1,nrow=3)
# ggsave(paste0(paste0(outPath,"/plots/",outputNamePrefix,"TADSinout_",
#                      padjVal,".pdf")),
#        width=9, height=11, paper="a4",plot=p,device="pdf")
#
#
#
# #separate anchors from inside tads
# tads_in<-reduce(GRanges(seqnames=seqnames(anchor1),IRanges(start=end(anchor1)+1,end=start(anchor2)-1)))
# tads_in<-resize(tads_in,width=width(tads_in)-20000,fix="center")
# anchors<-reduce(sort(c(anchor1,anchor2)))
# #anchors<-resize(anchors,width=width(anchors)+20000,fix="center")
# ol<-findOverlaps(anchors,tads_in)
# anchors<-anchors[-queryHits(ol)]
#
# plotList<-list()
# #grp=useContrasts[3]
# for (grp in useContrasts){
#   salmon<-readRDS(file=paste0(paste0(outPath,"/rds/",fileNamePrefix,
#                                      contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
#
#   salmon<-salmon[!is.na(salmon$chr),]
#   salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)
#
#   salmongr<-sort(salmongr)
#
#   ol<-findOverlaps(salmongr,tads_in,type="within")
#   insideTads<-salmongr[queryHits(ol)]
#
#   ol<-findOverlaps(salmongr,anchors)
#   atAnchors<-salmongr[queryHits(ol)]
#
#   insideTads$TADs<-"TAD"
#   atAnchors$TADs<-"Anchor"
#   df<-data.frame(c(insideTads,atAnchors))
#   df<-df%>%dplyr::group_by(seqnames,TADs)%>%dplyr::mutate(count=n())
#
#   plotList[[grp]]<-ggplot(df,aes(x=TADs,y=log2FoldChange,fill=TADs))+
#     geom_boxplot(notch=T,outlier.shape=NA,varwidth=T)+
#     facet_grid(.~seqnames) +ylim(c(-1,1))+
#     ggtitle(grp)
# }
# p<-gridExtra::marrangeGrob(plotList,ncol=1,nrow=3)
# ggsave(paste0(paste0(outPath,"/plots/",outputNamePrefix,"TADSvAnchors_",
#                      padjVal,".pdf")),
#        width=9, height=11, paper="a4",plot=p,device="pdf")
#
#
#
# ####################-
# ## new SIP loops -----
# ####################-
#
# # loops<-read.delim(paste0(outPath,"/otherData/10kbLoops.txt"),header=T)
# # gr1<-GRanges(seqnames=loops$chromosome1,ranges=IRanges(start=loops$x1,end=loops$x2-1))
# # gr2<-GRanges(seqnames=loops$chromosome2,ranges=IRanges(start=loops$y1,end=loops$y2-1))
# # mcols(gr1)<-loops[,c(8:15)]
# # mcols(gr2)<-loops[,c(8:15)]
# # gr1$loopNum<-paste0("loop",1:length(gr1))
# # gr2$loopNum<-paste0("loop",1:length(gr2))
# #
# # anchors<-c(gr1,gr2)
# # #anchors<-reduce(anchors,min.gapwidth=0L)
# # seqlevels(anchors)<-seqlevels(Celegans)[1:6]
# #
# # smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
# #                   useContrasts,"_lfc.bw")
# # names(smcRNAseq)<-useContrasts
# # for(grp in useContrasts){
# #   rnaSeq<-import.bw(smcRNAseq[[grp]])
# #   cov<-coverage(rnaSeq,weight="score")
# #   anchors<-binnedAverage(anchors,cov,grp)
# # }
# #
# #
# # pdf(paste0(outPath, "/plots/",outputNamePrefix,"RNAseq_loopsMetrics.pdf"),
# #            width=11,height=8,paper="a4r")
# # winSize=10000
# # metricsOI<-c("APScoreAvg","ProbabilityofEnrichment","RegAPScoreAvg","Avg_diffMaxNeihgboor_1", "Avg_diffMaxNeihgboor_2","avg","std","value")
# # plotList<-list()
# # for (colOI in metricsOI){
# #   df<-data.frame(anchors)
# #   colnames(df)<-c(colnames(df)[1:5],colnames(mcols(anchors)))
# #   df<-tidyr::pivot_longer(df,useContrasts,names_to="SMC",values_to="LFC")
# #   plotList[[colOI]]<-ggplot(df,aes_string(x=colOI,y="LFC",col="SMC")) + geom_point() +
# #     ylab("Average RNAseq score (10kb window)") + xlab(colOI) +
# #     ggtitle(paste0("Average RNAseq in ",winSize/1000,"kb window vs ",colOI))
# # }
# # p<-ggpubr::ggarrange(plotlist=plotList,ncol=2,nrow=2)
# # print(p)
# # dev.off()
#
#
#

#########################-
## insulation score ----
#########################-

insWin<-"500kb"
#insWin<-"250kb"
insulation<-read.delim(paste0(outPath,"/otherData/Insulation_10kb_",insWin,".csv"),
                       header=T,sep=" ")
#insulation<-read.delim(paste0(outPath,"/otherData/Insulation_10kb_250kb.csv"),
#                       header=T,sep=",")
ins.gr<-GRanges(seqnames=insulation$chrom, IRanges(start=insulation$start,
                                                   end=insulation$end),
                strand="*")
seqinfo(ins.gr)<-seqinfo(Celegans)
ins.gr$N2_L3<-insulation$N2_L3
ins.gr$N2_emb<-insulation$N2_Em
#export.bw(ins.gr,paste0(outPath,"/otherData/insulationScore_10kb_",insWin,".bw"))

listgr<-NULL
for (grp in useContrasts){
  #grp=useContrasts[1]
  salmon<-readRDS(file=paste0(outPath,"/rds/",fileNamePrefix,
                              contrastNames[[grp]],"_DESeq2_fullResults_p",
                              padjVal,".rds"))
  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)
  salmongr<-sort(salmongr)
  ol<-findOverlaps(salmongr,ins.gr)
  df<-data.frame(ol)
  df$N2_L3<-ins.gr$N2_L3[subjectHits(ol)]
  df$N2_emb<-ins.gr$N2_emb[subjectHits(ol)]
  df1<-df %>% group_by(queryHits) %>% summarise(l3=mean(N2_L3),emb=mean(N2_emb))
  salmongr$N2_L3<-NA
  salmongr$N2_L3[df1$queryHits]<-df1$l3
  salmongr$N2_emb<-NA
  salmongr$N2_emb[df1$queryHits]<-df1$emb
  salmongr$SMC<-grp
  listgr[[grp]]<-salmongr
}



insTbls<-do.call(rbind,lapply(listgr,data.frame))
row.names(insTbls)<-NULL
insTbls$XvA<-ifelse(insTbls$seqnames=="chrX","chrX","Autosomes")


p1<-ggplot(insTbls,aes(x=N2_L3,y=log2FoldChange)) +
  geom_point(size=1,color="#44444455") + ylim(c(-5,5))+
  facet_grid(SMC~XvA) +ggtitle(paste0(insWin," L3 insulation score Vs LFC"))
p1
ggsave(paste0(outPath, "/plots/",outputNamePrefix,"LFCvsL3InsulationScore.png"),
       p1, device="png",width=19,height=29,units="cm")


p2<-ggplot(insTbls,aes(x=N2_emb,y=log2FoldChange)) +
  geom_point(size=1,color="#44444455") + ylim(c(-5,5))+
  facet_grid(SMC~XvA) +ggtitle(paste0(insWin," Emb insulation score Vs LFC"))
p2
ggsave(paste0(outPath, "/plots/",outputNamePrefix,"LFCvsEmbInsulationScore.png"),
       p2, device="png",width=19,height=29,units="cm")




# ##########################-
# ########### moustache loops------
# ##########################-
#
# ### new Moustache loops
# #mustacheBatch="PMW366"
# mustacheBatch="PMW382"
#
# loops<-import(paste0(outPath,"/otherData/",mustacheBatch,"_2k_mustache_filtered.bedpe"),format="bedpe")
# grl<-zipup(loops)
# anchor1<-do.call(c,lapply(grl,"[",1))
# anchor2<-do.call(c,lapply(grl,"[",2))
# mcols(anchor1)<-mcols(loops)
# mcols(anchor2)<-mcols(loops)
#
# anchor1$loopNum<-paste0("loop",1:length(anchor1))
# anchor2$loopNum<-paste0("loop",1:length(anchor2))
#
# # anchors<-c(anchor1,anchor2)
# # #anchors<-reduce(anchors,min.gapwidth=0L)
# # seqlevels(anchors)<-seqlevels(Celegans)[1:6]
#
# #make TADs
# tads<-GRanges(seqnames=seqnames(anchor1),IRanges(start=start(anchor1),end=end(anchor2)))
# head(tads)
# sort(width(tads))
# tads<-reduce(tads)
# sort(width(tads))
#
# # find regions not in tads
# notads<-gaps(tads)
# sort(width(notads))
#
#
# plotList<-list()
# #grp=useContrasts[3]
# for (grp in useContrasts){
#   salmon<-readRDS(file=paste0(paste0(outPath,"/rds/",fileNamePrefix,
#                                      contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
#
#   salmon<-salmon[!is.na(salmon$chr),]
#   salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)
#
#   salmongr<-sort(salmongr)
#
#   ol<-findOverlaps(salmongr,tads,type="within")
#   genesInTads<-salmongr[queryHits(ol)]
#
#   ol<-findOverlaps(salmongr,notads)
#   genesNotTads<-salmongr[queryHits(ol)]
#
#   genesInTads$Loops<-"inside"
#   genesNotTads$Loops<-"outside"
#   df<-data.frame(c(genesInTads,genesNotTads))
#   df<-df%>%dplyr::group_by(seqnames,Loops)%>%dplyr::mutate(count=n())
#
#   plotList[[grp]]<-ggplot(df,aes(x=Loops,y=log2FoldChange,fill=Loops))+
#     geom_boxplot(notch=T,outlier.shape=NA,varwidth=T)+
#     facet_grid(.~seqnames) +ylim(c(-1,1))+
#     ggtitle(grp)
# }
# p<-gridExtra::marrangeGrob(plotList,ncol=1,nrow=3)
# ggsave(paste0(paste0(outPath,"/plots/",outputNamePrefix,"LoopsInOut_",mustacheBatch,"Mustache_",
#                      padjVal,".pdf")),
#        width=9, height=11, paper="a4",plot=p,device="pdf")
#
#
#
# #separate anchors from inside tads
# tads_in<-GRanges(seqnames=seqnames(anchor1),IRanges(start=end(anchor1)+1,end=start(anchor2)-1))
# # tads_in<-reduce(tads_in)
# tads_in<-resize(tads_in,width=width(tads_in)-20000,fix="center")
# anchors<-sort(c(anchor1,anchor2))
# anchors<-resize(anchors,width=20000,fix="center")
# #anchors<-reduce(anchors)
# #ol<-findOverlaps(anchors,tads_in)
# #anchors<-anchors[-queryHits(ol)]
#
# width(anchors)
# dataList<-list()
# plotList<-list()
# #grp=useContrasts[3]
# for (grp in useContrasts){
#   salmon<-readRDS(file=paste0(paste0(outPath,"/rds/",fileNamePrefix,
#                                      contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
#
#   salmon<-salmon[!is.na(salmon$chr),]
#   salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)
#
#   salmongr<-sort(salmongr)
#
#   ol<-findOverlaps(salmongr,tads_in,type="within")
#   insideTads<-salmongr[queryHits(ol)]
#
#   ol<-findOverlaps(salmongr,anchors)
#   atAnchors<-salmongr[queryHits(ol)]
#
#   insideTads$Loops<-"inLoop"
#   atAnchors$Loops<-"Anchor"
#   df<-data.frame(c(insideTads,atAnchors))
#   df<-df%>%dplyr::group_by(seqnames,Loops)%>%dplyr::mutate(count=n())
#   df$SMC<-grp
#
#   dataList[[grp]]<-df
#   plotList[[grp]]<-ggplot(df,aes(x=Loops,y=log2FoldChange,fill=Loops))+
#     geom_boxplot(notch=T,outlier.shape=NA,varwidth=T)+
#     facet_grid(.~seqnames) +ylim(c(-1,1))+
#     ggtitle(grp)
# }
# p<-gridExtra::marrangeGrob(plotList,ncol=1,nrow=3)
#
# ggsave(paste0(paste0(outPath,"/plots/",outputNamePrefix,"LoopsvAnchors_",mustacheBatch,"-Mostache_",
#                      padjVal,".pdf")),
#        width=9, height=11, paper="a4",plot=p,device="pdf")
#
# ## focus on chrX loops
# dataTbl<-do.call(rbind,dataList)
# xchr<-dataTbl[dataTbl$seqnames=="chrX",]
#
# xchr$SMC<-factor(xchr$SMC,levels=useContrasts)
# p1<-ggplot(xchr,aes(x=Loops,y=log2FoldChange,fill=Loops))+
#   geom_boxplot(notch=T,outlier.shape=NA,varwidth=T)+
#   facet_grid(~SMC) +ylim(c(-1,1))+
#   ggtitle(paste0("LFC at ",mustacheBatch," anchors Vs inside loops in chrX")) +
#   geom_hline(yintercept=0,linetype="dotted",color="grey20") +
#   theme(axis.text.x=element_text(angle=45,hjust=1))+
#   xlab(label=element_blank())
#
# xchr$SMC<-factor(xchr$SMC,levels=useContrasts)
# xchr$measure="Expression"
# p2<-ggplot(xchr,aes(x=Loops,y=log2(baseMean),fill=Loops))+
#   geom_boxplot(notch=T,outlier.shape=NA,varwidth=T) +
#   facet_wrap(.~measure) + ggtitle("Base mean counts") +
#   theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1)) +
#   xlab(label=element_blank())
#
# p<-ggarrange(p2,p1,ncol=2,widths=c(1.3,8.7))
# ggsave(paste0(paste0(outPath,"/plots/",outputNamePrefix,"LoopsvAnchorsXchr_",mustacheBatch,"-Mostache_",
#                      padjVal,".pdf")),
#        width=11, height=5, paper="a4r",plot=p,device="pdf")
#
#
#
# ################-
# ## compare anchors
# ################-
#
# loops366<-import(paste0(outPath,"/otherData/PMW366_2k_mustache_filtered.bedpe"),format="bedpe")
#
# loops382<-import(paste0(outPath,"/otherData/PMW382_2k_mustache_filtered.bedpe"),format="bedpe")
#
#

##########################-
## Manual clicked loops------
##########################-

### new clicked loops
#clickedBatch="366"
clickedBatch="382"

ceTiles<-tileGenome(seqlengths(Celegans),tilewidth=10000,cut.last.tile.in.chrom = T)

loops<-import(paste0(outPath,"/otherData/Clicked_loops_",clickedBatch,"_merge.bedpe"),format="bedpe")
grl<-zipup(loops)
anchor1<-do.call(c,lapply(grl,"[",1))
anchor2<-do.call(c,lapply(grl,"[",2))
mcols(anchor1)<-mcols(loops)
mcols(anchor2)<-mcols(loops)

anchor1$loopNum<-paste0("loop",1:length(anchor1))
anchor2$loopNum<-paste0("loop",1:length(anchor2))

#separate anchors from inside tads
tads_in<-GRanges(seqnames=seqnames(anchor1),IRanges(start=end(anchor1)+1,end=start(anchor2)-1))
# tads_in<-reduce(tads_in)
tads_in<-resize(tads_in,width=width(tads_in)-20000,fix="center")

ol<-findOverlaps(ceTiles,tads_in)
tenkbInTads<-ceTiles[unique(queryHits(ol))]

anchors<-sort(c(anchor1,anchor2))
anchors<-resize(anchors,width=10000,fix="center")
reduce(anchors)

ol<-findOverlaps(tenkbInTads,anchors)
tenkbInTads<-tenkbInTads[-queryHits(ol)]

width(tenkbInTads)
width(anchors)
dataList<-list()
plotList<-list()
#grp=useContrasts[3]
for (grp in useContrasts){
  salmon<-readRDS(file=paste0(paste0(outPath,"/rds/",fileNamePrefix,
                                     contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))

  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

  salmongr<-sort(salmongr)

  ol<-findOverlaps(salmongr,tenkbInTads)
  insideTads<-salmongr[queryHits(ol)]

  ol<-findOverlaps(salmongr,anchors)
  atAnchors<-salmongr[queryHits(ol)]

  insideTads$Loops<-"inLoop"
  atAnchors$Loops<-"Anchor"
  df<-data.frame(c(insideTads,atAnchors))
  df<-df%>%dplyr::group_by(seqnames,Loops)%>%dplyr::mutate(count=n())
  df$SMC<-grp

  dataList[[grp]]<-df
  plotList[[grp]]<-ggplot(df,aes(x=Loops,y=log2FoldChange,fill=Loops))+
    geom_boxplot(notch=T,outlier.shape=NA,varwidth=T)+
    facet_grid(.~seqnames) +ylim(c(-1,1))+
    ggtitle(grp)
}
p<-gridExtra::marrangeGrob(plotList,ncol=1,nrow=3)

ggsave(paste0(paste0(outPath,"/plots/",outputNamePrefix,"LoopsvAnchors_",clickedBatch,"-Clicked_",
                     padjVal,".pdf")),
       width=9, height=11, paper="a4",plot=p,device="pdf")

## focus on chrX loops
dataTbl<-do.call(rbind,dataList)
xchr<-dataTbl[dataTbl$seqnames=="chrX",]

xchr$SMC<-factor(xchr$SMC,levels=useContrasts)
p1<-ggplot(xchr,aes(x=Loops,y=log2FoldChange,fill=Loops))+
  geom_boxplot(notch=T,outlier.shape=NA,varwidth=T)+
  facet_grid(~SMC) +ylim(c(-1,1))+
  ggtitle(paste0("LFC at ",clickedBatch," anchors Vs inside loops in chrX")) +
  geom_hline(yintercept=0,linetype="dotted",color="grey20") +
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  xlab(label=element_blank())

xchr$SMC<-factor(xchr$SMC,levels=useContrasts)
xchr$measure="Expression"
p2<-ggplot(xchr,aes(x=Loops,y=log2(baseMean),fill=Loops))+
  geom_boxplot(notch=T,outlier.shape=NA,varwidth=T) +
  facet_wrap(.~measure) + ggtitle("Base mean counts") +
  theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1)) +
  xlab(label=element_blank())

p<-ggarrange(p2,p1,ncol=2,widths=c(1.3,8.7))
ggsave(paste0(paste0(outPath,"/plots/",outputNamePrefix,"LoopsvAnchorsXchr_",clickedBatch,"-Mostache_",
                     padjVal,".pdf")),
       width=11, height=5, paper="a4r",plot=p,device="pdf")



###################-
## binned signal around anchors
###################-

anchordf<-data.frame(source=c("clicked366","clicked382","eigen382"),
                     file=c(paste0(outPath,"/otherData/all_anchors_loops_",
                                   c("366","382"),"_full_size_correct.bed"),
                            paste0(outPath,"/otherData/382_X.eigs_cis.vecs_37peaks_p0.65_correct.bed")))
### LFC --------
contrastsToUse<-useContrasts[c(1,3,4,6,7,8,9)]
smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
                  contrastsToUse,"_lfc.bw")
names(smcRNAseq)<-contrastsToUse

for(anch in 1:nrow(anchordf)){
  anchors<-import(anchordf$file[anch])
  anchorSource<-anchordf$source[anch]

  anchors<-anchors[seqnames(anchors)=="chrX"]
  seqlevels(anchors)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
  seqinfo(anchors)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)

  winSize=10000
  anchors<-resize(anchors,width=winSize,fix="center")
  seqlevels(anchors)<-seqlevels(Celegans)[1:6]

  xanchors<-anchors[seqnames(anchors)=="chrX"]

  p<-avrSignalBins(xanchors, bwFiles=smcRNAseq, winSize=winSize,numWins=100000/winSize)
  p<-p+ggtitle(paste0("Average RNAseq LFC around chrX ",anchorSource," loop anchors"))
  print(p)
  ggsave(paste0(outPath, "/plots/",outputNamePrefix,"binnedRNAseq_",anchorSource,"anchors_win",
                winSize/1000,"kb.pdf"), p,
         device="pdf",width=19,height=29,units="cm")
}

### LFC unfiltered------
if(dir.exists(paste0(outPath,"/tracks/p0.05_lfc0.5"))){
  localFileNamePrefix="p0.05_lfc0.5/no775B3_"
  contrastsToUse<-useContrasts[c(3,1,4,6,7,8,9)]
  smcRNAseq<-paste0(outPath,"/tracks/",localFileNamePrefix,
                    contrastsToUse,"_wt_lfc.bw")
  names(smcRNAseq)<-contrastsToUse

  for(anch in 1:nrow(anchordf)){
    anchors<-import(anchordf$file[anch])
    anchorSource<-anchordf$source[anch]

    anchors<-anchors[seqnames(anchors)=="chrX"]
    seqlevels(anchors)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
    seqinfo(anchors)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)

    winSize=10000
    anchors<-resize(anchors,width=winSize,fix="center")
    seqlevels(anchors)<-seqlevels(Celegans)[1:6]

    xanchors<-anchors[seqnames(anchors)=="chrX"]

    p<-avrSignalBins(xanchors, bwFiles=smcRNAseq, winSize=winSize,numWins=100000/winSize)
    p<-p+ggtitle(paste0("Average RNAseq LFC around chrX ",anchorSource," loop anchors"))
    print(p)
    ggsave(paste0(outPath, "/plots/",outputNamePrefix,"binnedRNAseq_",anchorSource,"anchors_win",
                  winSize/1000,"kb_LFCunfiltered.pdf"), p,
           device="pdf",width=19,height=29,units="cm")
  }
}


### TPM -------
anchordf<-data.frame(source=c("clicked366","clicked382","eigen382"),
                     file=c(paste0(outPath,"/otherData/all_anchors_loops_",
                                   c("366","382"),"_full_size_correct.bed"),
                            paste0(outPath,"/otherData/382_X.eigs_cis.vecs_37peaks_p0.65_correct.bed")))

contrastsToUse<-c("366","382","775","784","828","822")
smcRNAseq<-paste0(outPath,"/tracks/PMW",contrastsToUse,"_TPM_avr.bw")
names(smcRNAseq)<-paste0("PMW",contrastsToUse)

for(anch in 1:nrow(anchordf)){
  anchors<-import(anchordf$file[anch])
  anchorSource<-anchordf$source[anch]

  anchors<-anchors[seqnames(anchors)=="chrX"]
  seqlevels(anchors)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
  seqinfo(anchors)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)

  winSize=10000
  anchors<-resize(anchors,width=winSize,fix="center")
  seqlevels(anchors)<-seqlevels(Celegans)

  xanchors<-anchors[seqnames(anchors)=="chrX"]

  p<-avrSignalBins(xanchors, bwFiles=smcRNAseq, winSize=winSize,
                   numWins=100000/winSize, logScore=T)
  p<-p+ggtitle(paste0("Average RNAseq TPM around chrX ",anchorSource," loop anchors"))
  print(p)
  ggsave(paste0(outPath, "/plots/",outputNamePrefix,"RNAseqTPM_",anchorSource,"anchors_win",
                winSize/1000,"kb.pdf"), p,
         device="pdf",width=19,height=29,units="cm")
}



###############-
## seqplots around anchors-------
###############-

flankSize<-100000

anchordf<-data.frame(source=c("366","382","eigen382"),
                     file=c(paste0(outPath,"/otherData/all_anchors_loops_",
                                   c("366","382"),"_full_size_correct.bed"),
                            paste0(outPath,"/otherData/382_X.eigs_cis.vecs_37peaks_p0.65_correct.bed")))

### LFC ------
contrastsToUse<-useContrasts[c(3,1,4,6,7,8,9)]
smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
                  contrastsToUse,"_lfc.bw")
names(smcRNAseq)<-contrastsToUse


for(anch in 1:nrow(anchordf)){
  anchors<-import(anchordf$file[anch])
  anchorSource<-anchordf$source[anch]

  anchors<-anchors[seqnames(anchors)=="chrX"]
  seqlevels(anchors)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
  seqinfo(anchors)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)

  seqlevels(anchors)<-seqlevels(Celegans)[1:6]

  xanchors<-anchors[seqnames(anchors)=="chrX"]

  pdf(file=paste0(outPath,"/plots/",outputNamePrefix,anchorSource,"anchors_flank",
                  flankSize/1000,"kb_LFC.pdf"), width=11,
      height=9, paper="a4r")


  p<-getPlotSetArray(tracks=c(smcRNAseq),
                     features=c(xanchors),
                     refgenome="ce11", bin=flankSize/100, xmin=flankSize,
                     xmax=flankSize, type="af", rm0=T, ignore_strand=T,
                     xanchored=median(width(xanchors)))



  dd<-plotHeatmap(p,plotz=F)
  heatmapQuantiles<-sapply(dd$HLST,quantile,c(0.01,0.99),na.rm=T)
  roworder<-rev(order(lapply(dd$HLST,rowSums,na.rm=T)$X1))
  minVal<-min(heatmapQuantiles[1,])
  maxVal<-max(heatmapQuantiles[2,])
  #layout(matrix(c(1), nrow = 2, ncol = 1, byrow = TRUE))
  plotAverage(p,ylim=c(-0.5,1),main="ChrX anchors",
              error.estimates=ifelse(length(useContrasts>3),F,T))
  plotHeatmap(p,main="ChrX anchors", plotScale="linear", sortrows=T,
              clusters=1L,autoscale=F,zmin=minVal, zmax=maxVal, ln.v=F,
              indi=F, sort_mids=T,sort_by=c(T,rep(F,length(smcRNAseq)-1)),
              clspace=c("#00008B", "#FFFFE0","#8B0000"))
  dev.off()
}



### LFC unfiltered------
if(dir.exists(paste0(outPath,"/tracks/p0.05_lfc0.5"))){
  localFileNamePrefix="p0.05_lfc0.5/no775B3_"
  contrastsToUse<-useContrasts[c(3,1,4,6,7,8,9)]
  smcRNAseq<-paste0(outPath,"/tracks/",localFileNamePrefix,
                    contrastsToUse,"_wt_lfc.bw")
  names(smcRNAseq)<-contrastsToUse

  for(anch in 1:nrow(anchordf)){
    anchors<-import(anchordf$file[anch])
    anchorSource<-anchordf$source[anch]

    anchors<-anchors[seqnames(anchors)=="chrX"]
    seqlevels(anchors)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
    seqinfo(anchors)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)

    seqlevels(anchors)<-seqlevels(Celegans)#[1:6]

    xanchors<-anchors[seqnames(anchors)=="chrX"]

    pdf(file=paste0(outPath,"/plots/",outputNamePrefix,anchorSource,"anchors_flank",
                    flankSize/1000,"kb_LFCunfilt.pdf"), width=11,
        height=9, paper="a4r")

    p<-getPlotSetArray(tracks=c(smcRNAseq),
                       features=c(xanchors),
                       refgenome="ce11", bin=flankSize/100, xmin=flankSize,
                       xmax=flankSize, type="af", rm0=T, ignore_strand=T,
                       xanchored=median(width(xanchors)))

    dd<-plotHeatmap(p,plotz=F)
    heatmapQuantiles<-sapply(dd$HLST,quantile,c(0.01,0.99),na.rm=T)
    roworder<-rev(order(lapply(dd$HLST,rowSums,na.rm=T)$X1))
    minVal<-min(heatmapQuantiles[1,])
    maxVal<-max(heatmapQuantiles[2,])
    #layout(matrix(c(1), nrow = 2, ncol = 1, byrow = TRUE))
    plotAverage(p,ylim=c(-0.5,1),main="ChrX anchors",
                error.estimates=ifelse(length(useContrasts>3),F,T))
    plotHeatmap(p,main="ChrX anchors", plotScale="linear", sortrows=T,
                clusters=1L,autoscale=F,zmin=minVal, zmax=maxVal,
                ln.v=F,
                indi=F, sort_mids=T,sort_by=c(T,rep(F,length(smcRNAseq)-1)),
                clspace=c("#00008B", "#FFFFE0","#8B0000"))
    dev.off()
  }
}


#### TPM ------
anchordf<-data.frame(source=c("clicked366","clicked382","eigen382"),
                     file=c(paste0(outPath,"/otherData/all_anchors_loops_",
                                   c("366","382"),"_full_size_correct.bed"),
                            paste0(outPath,"/otherData/382_X.eigs_cis.vecs_37peaks_p0.65_correct.bed")))

contrastsToUse<-c("366","382","775","784","828","822")
smcRNAseq<-paste0(outPath,"/tracks/PMW",contrastsToUse,"_TPM_avr.bw")
names(smcRNAseq)<-paste0("PMW",contrastsToUse)

for(anch in 1:nrow(anchordf)){
  anchors<-import(anchordf$file[anch])
  anchorSource<-anchordf$source[anch]

  anchors<-anchors[seqnames(anchors)=="chrX"]
  seqlevels(anchors)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
  seqinfo(anchors)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)

  seqlevels(anchors)<-seqlevels(Celegans)[1:6]

  xanchors<-anchors[seqnames(anchors)=="chrX"]

  pdf(file=paste0(outPath,"/plots/",outputNamePrefix,anchorSource,"anchors_flank",
                  flankSize/1000,"kb_TPM.pdf"), width=11,
      height=9, paper="a4r")


  p<-getPlotSetArray(tracks=c(smcRNAseq),
                     features=c(xanchors),
                     refgenome="ce11", bin=flankSize/10, xmin=flankSize,
                     xmax=flankSize, type="af", rm0=T, ignore_strand=T,
                     xanchored=median(width(xanchors)))

  dd<-plotHeatmap(p,plotz=F)
  heatmapQuantiles<-sapply(dd$HLST,quantile,c(0.01,0.99),na.rm=T)
  roworder<-rev(order(lapply(dd$HLST,rowSums,na.rm=T)$X1))
  minVal<-min(heatmapQuantiles[1,])
  maxVal<-max(heatmapQuantiles[2,])
  #layout(matrix(c(1), nrow = 2, ncol = 1, byrow = TRUE))
  plotAverage(p,main="ChrX anchors", plotScale="log2",
              error.estimates=ifelse(length(useContrasts>3),F,T))
  plotHeatmap(p,main="ChrX anchors", plotScale="log2", sortrows=T,
              clusters=1L,autoscale=F,
              zmin=minVal, zmax=maxVal,
              ln.v=F,
              indi=F, sort_mids=T,sort_by=c(T,rep(F,length(smcRNAseq)-1)),
              clspace=c("#00008B", "#FFFFE0","#8B0000"))
  dev.off()
}





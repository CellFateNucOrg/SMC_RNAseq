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
pca1<-import.bw(paste0(outPath,"/otherData/N2_merge_2000.E1.vecs.bw"))
pca2<-import.bw(paste0(outPath,"/otherData/N2_merge_2000.E2.vecs.bw"))

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

#' Collect counts of significantly changed genes per chromosome and per compartment
#'
#' @param listgr List of GRanges for different RNAseq
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


plotCountsPerChrPerCompartment<-function(dfl,namePCAcol,namePCA){
  ymax=max(dfl$n)
  p<-ggplot(dfl,aes(x=seqnames,y=n,group=get(namePCAcol))) +
    geom_bar(stat="identity", position=position_dodge(), aes(fill=get(namePCAcol))) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
    xlab("chr")+ylab("Number of genes") +
    theme(legend.title = element_text(size=10)) +
    ggtitle(paste0("Significantly changed genes per chromosome by ",namePCA," compartment")) + labs(fill=namePCA)
  return(p)
}



plotFractionPerChrPerCompartment<-function(dfl,namePCAcol,namePCA){
  p<-ggplot(dfl,aes(x=seqnames,y=Frac,group=get(namePCAcol))) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=get(namePCAcol))) +
  facet_grid(cols=vars(SMC)) +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle(paste0("Fraction changed genes per chromosome by ",namePCA," compartment")) + labs(fill=namePCA)
  return(p)
}


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


plotUpDownByCompartment<-function(dfl,namePCAcol,namePCA,compartment){
  yminmax=c(0,max(dfl$n))
  p<-ggplot(dfl[dfl[,namePCAcol]==compartment,],aes(x=seqnames,y=n,group=expression)) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=expression)) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
    xlab("chr")+ylab("Number of genes") +
    ggtitle(paste0("Up/down regulated in ",namePCA," ",compartment,
                   " compartment"))
  return(p)
}

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

localLFC=0
#localLFC=lfcVal

# upregulated
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=localLFC,direction="gt")

sigList<-lapply(sigList, "[", ,c("N2.E1_compartment","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList<-sigList[!is.na(sigList$N2.E1_compartment),]
sigTbl<-sigList
sigTbl$updown<-"up"


# downregulated
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal, lfc= -localLFC, direction="lt")

sigList<-lapply(sigList, "[", ,c("N2.E1_compartment","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList<-sigList[!is.na(sigList$N2.E1_compartment),]
sigList$updown<-"down"
sigTbl<-rbind(sigTbl,sigList)
sigTbl$N2.E1_compartment<-as.factor(sigTbl$N2.E1_compartment)
sigTbl$updown<-factor(sigTbl$updown, levels=c("up","down"))

yminmax=c(0,max(abs(sigTbl$log2FoldChange)))
p2<-ggplot(sigTbl,aes(x=N2.E1_compartment,y=abs(log2FoldChange),col=updown,fill=updown)) +
  geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2),outlier.size=0.4,
               outlier.color="grey50") +
  facet_grid(cols=vars(SMC)) +
  ggtitle("Significantly changed genes by N2 E1 compartment") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide="none")

yminmax=c(0,median(abs(sigTbl$log2FoldChange))+quantile(abs(sigTbl$log2FoldChange))[4]*2)
p3<-ggplot(sigTbl,aes(x=N2.E1_compartment,y=abs(log2FoldChange),col=updown,fill=updown)) +
  geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2), outlier.shape=NA,
               outlier.color="grey50") +
  facet_grid(cols=vars(SMC)) +
  ggtitle("Significantly changed genes by N2 E1 compartment") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide="none")

# test significance of LFC
summary(aov(abs(log2FoldChange)~updown+N2.E1_compartment,data=sigTbl))

p<-ggpubr::ggarrange(p2,p3,ncol=1,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                "ABcomp_N2-E1_LFC_padj",
                          padjVal,"_lfc", localLFC,".pdf"),
                plot=p, device="pdf",width=29,height=19,units="cm")

### second N2 eigenvector------

# upregulated
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=localLFC,direction="gt")

sigList<-lapply(sigList, "[", ,c("N2.E2_compartment","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList<-sigList[!is.na(sigList$N2.E2_compartment),]
sigTbl<-sigList
sigTbl$updown<-"up"


# downregulated
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal, lfc= -localLFC, direction="lt")

sigList<-lapply(sigList, "[", ,c("N2.E2_compartment","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList<-sigList[!is.na(sigList$N2.E2_compartment),]
sigList$updown<-"down"
sigTbl<-rbind(sigTbl,sigList)
sigTbl$N2.E2_compartment<-as.factor(sigTbl$N2.E2_compartment)
sigTbl$updown<-factor(sigTbl$updown, levels=c("up","down"))

yminmax=c(0,max(abs(sigTbl$log2FoldChange)))
p2<-ggplot(sigTbl,aes(x=N2.E2_compartment,y=abs(log2FoldChange),col=updown,fill=updown)) +
  geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2),outlier.size=0.4,
               outlier.color="grey50") +
  facet_grid(cols=vars(SMC)) +
  ggtitle("Significantly changed genes by N2 E2 compartment") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide="none")

yminmax=c(0,median(abs(sigTbl$log2FoldChange))+quantile(abs(sigTbl$log2FoldChange))[4]*2)
p3<-ggplot(sigTbl,aes(x=N2.E2_compartment,y=abs(log2FoldChange),col=updown,fill=updown)) +
  geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2), outlier.shape=NA,
               outlier.color="grey50") +
  facet_grid(cols=vars(SMC)) +
  ggtitle("Significantly changed genes by N2 E2 compartment") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide="none")

# test significance of LFC
summary(aov(abs(log2FoldChange)~updown+N2.E2_compartment,data=sigTbl))

p<-ggpubr::ggarrange(p2,p3,ncol=1,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                "ABcomp_N2-E2_LFC_padj",
                                padjVal,"_lfc", localLFC,".pdf"),
                plot=p, device="pdf",width=29,height=19,units="cm")






####-
## AB compartment by chromosome 366-----
####-

####
## 366 compartments
####
pca1<-import.bw(paste0(outPath,"/otherData/366_merge_2000.E1.vecs.bw"))
pca2<-import.bw(paste0(outPath,"/otherData/366_merge_2000.E2.vecs.bw"))


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

######### 366 first Eigenvector ------
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


########## 366 second eigenvector ------
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
#### first 366 eigenvector------
#localLFC=0
localLFC=lfcVal

# upregulated
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=localLFC,direction="gt")

sigList<-lapply(sigList, "[", ,c("PMW366.E1_compartment","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList<-sigList[!is.na(sigList$PMW366.E1_compartment),]
sigTbl<-sigList
sigTbl$updown<-"up"


# downregulated
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal, lfc= -localLFC, direction="lt")

sigList<-lapply(sigList, "[", ,c("PMW366.E1_compartment","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList<-sigList[!is.na(sigList$PMW366.E1_compartment),]
sigList$updown<-"down"
sigTbl<-rbind(sigTbl,sigList)
sigTbl$PMW366.E1_compartment<-as.factor(sigTbl$PMW366.E1_compartment)
sigTbl$updown<-factor(sigTbl$updown, levels=c("up","down"))

yminmax=c(0,max(abs(sigTbl$log2FoldChange)))
p2<-ggplot(sigTbl,aes(x=PMW366.E1_compartment,y=abs(log2FoldChange),col=updown,fill=updown)) +
  geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2),outlier.size=0.4,
               outlier.color="grey50") +
  facet_grid(cols=vars(SMC)) +
  ggtitle("Significantly changed genes by 366 E1 compartment") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
  theme(legend.title = element_text(size=10)) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide="none")

yminmax=c(0,median(abs(sigTbl$log2FoldChange))+quantile(abs(sigTbl$log2FoldChange))[4]*2)
p3<-ggplot(sigTbl,aes(x=PMW366.E1_compartment,y=abs(log2FoldChange),col=updown,fill=updown)) +
  geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2), outlier.shape=NA,
               outlier.color="grey50") +
  facet_grid(cols=vars(SMC)) +
  ggtitle("Significantly changed genes by 366 E1 compartment") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
  theme(legend.title = element_text(size=10)) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide="none")

# test significance of LFC
summary(aov(abs(log2FoldChange)~updown+PMW366.E1_compartment,data=sigTbl))

p<-ggpubr::ggarrange(p2,p3,ncol=1,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                "ABcomp_PMW366-E1_LFC_padj",
                                padjVal,"_lfc", localLFC,".pdf"),
                plot=p, device="pdf",width=29,height=19,units="cm")

### second 366 eigenvector------
# genes that change significantly
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=localLFC,direction="both")
#sigList<-lapply(listgr, as.data.frame)

sigList<-lapply(listgr,as.data.frame)

sigList<-lapply(sigList, "[", ,c("PMW366.E2_compartment","log2FoldChange"))
# #sigList$SMC<-NA
# for(g in names(sigList)){ sigList[[g]]$SMC<-g }
# sigList<-do.call(rbind,sigList)
# sigList<-sigList[!is.na(sigList$compartment),]
# sigList$compartment<-as.factor(sigList$compartment)
#
# yminmax=max(abs(min(sigList$log2FoldChange)),max(sigList$log2FoldChange))
# yminmax<-c(-yminmax,yminmax)
# p1<-ggplot(sigList,aes(x=compartment,y=log2FoldChange,fill=compartment)) +
#   geom_violin() + facet_grid(cols=vars(SMC)) +
#   ylim(yminmax) +
#   ggtitle("Significantly changed genes by N2 compartment") +
#   theme_minimal() + scale_fill_grey(start=0.8,end=0.3)


# upregulated
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=localLFC,direction="gt")

sigList<-lapply(sigList, "[", ,c("PMW366.E2_compartment","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList<-sigList[!is.na(sigList$PMW366.E2_compartment),]
sigTbl<-sigList
sigTbl$updown<-"up"


# downregulated
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal, lfc= -localLFC, direction="lt")

sigList<-lapply(sigList, "[", ,c("PMW366.E2_compartment","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList<-sigList[!is.na(sigList$PMW366.E2_compartment),]
sigList$updown<-"down"
sigTbl<-rbind(sigTbl,sigList)
sigTbl$PMW366.E2_compartment<-as.factor(sigTbl$PMW366.E2_compartment)
sigTbl$updown<-factor(sigTbl$updown, levels=c("up","down"))

yminmax=c(0,max(abs(sigTbl$log2FoldChange)))
p2<-ggplot(sigTbl,aes(x=PMW366.E2_compartment,y=abs(log2FoldChange),col=updown,fill=updown)) +
  geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2),outlier.size=0.4,
               outlier.color="grey50") +
  facet_grid(cols=vars(SMC)) +
  ggtitle("Significantly changed genes by 366 E2 compartment") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
  theme(legend.title = element_text(size=10)) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide="none")

yminmax=c(0,median(abs(sigTbl$log2FoldChange))+quantile(abs(sigTbl$log2FoldChange))[4]*2)
p3<-ggplot(sigTbl,aes(x=PMW366.E2_compartment,y=abs(log2FoldChange),col=updown,fill=updown)) +
  geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2), outlier.shape=NA,
               outlier.color="grey50") +
  facet_grid(cols=vars(SMC)) +
  ggtitle("Significantly changed genes by 366 E2 compartment") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
  theme(legend.title = element_text(size=10)) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide="none")

# test significance of LFC
summary(aov(abs(log2FoldChange)~updown+PMW366.E2_compartment,data=sigTbl))

p<-ggpubr::ggarrange(p2,p3,ncol=1,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                "ABcomp_PMW366-E2_LFC_padj",
                                padjVal,"_lfc", localLFC,".pdf"),
                plot=p, device="pdf",width=29,height=19,units="cm")


#############-
# Both eigen vectors combined------
#############-

localLFC=0
#localLFC=lfcVal

# unfiltered
sigList<-lapply(lapply(listgr,as.data.frame), "[", ,c("wormbaseID","PMW366.E1_compartment","PMW366.E2_compartment","log2FoldChange","padj"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList$PMW366.E1.E2<-paste0(sigList$PMW366.E1_compartment,"1.",sigList$PMW366.E2_compartment,"2")
sigList<-sigList[!is.na(sigList$padj),]
sigList$PMW366.E1.E2<-factor(sigList$PMW366.E1.E2)
totalPerCat<-sigList %>% group_by(SMC,PMW366.E1.E2) %>% summarise(count=n())

ggplot(totalPerCat,aes(x=PMW366.E1.E2,y=count,group=SMC))+geom_bar(stat="identity")+
  facet_wrap(vars(SMC))

#yminmax=c(0,median(abs(sigList$log2FoldChange))+quantile(abs(sigList$log2FoldChange))[4]*2)
yminmax=c(0,0.6)
give.n <- function(x){
  return(c(y = 0.5, label = length(x)))
  # experiment with the multiplier to find the perfect position
}
ggplot(sigList,aes(x=PMW366.E1.E2,y=abs(log2FoldChange))) +
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

sigList<-lapply(sigList, "[", ,c("PMW366.E1_compartment","PMW366.E2_compartment","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList$PMW366.E1.E2<-paste0(sigList$PMW366.E1_compartment,"1.",sigList$PMW366.E2_compartment,"2")

sigList<-sigList[!is.na(sigList$PMW366.E1.E2),]
sigTbl<-sigList
sigTbl$updown<-"up"


# downregulated
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal, lfc= -localLFC, direction="lt")

sigList<-lapply(sigList, "[", ,c("PMW366.E1_compartment","PMW366.E2_compartment","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList$PMW366.E1.E2<-paste0(sigList$PMW366.E1_compartment,"1.",sigList$PMW366.E2_compartment,"2")
sigList<-sigList[!is.na(sigList$PMW366.E1.E2),]
sigList$updown<-"down"
sigTbl<-rbind(sigTbl,sigList)
sigTbl$PMW366.E1.E2<-as.factor(sigTbl$PMW366.E1.E2)
sigTbl$updown<-factor(sigTbl$updown, levels=c("up","down"))

yminmax=c(0,median(abs(sigTbl$log2FoldChange))+quantile(abs(sigTbl$log2FoldChange))[4]*2)
give.n <- function(x){
  return(c(y = 0, label = length(x)))
  # experiment with the multiplier to find the perfect position
}

p2<-ggplot(sigTbl,aes(x=PMW366.E1.E2,y=abs(log2FoldChange),col=updown,fill=updown)) +
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
summary(aov(abs(log2FoldChange)~updown+PMW366.E1.E2,data=sigTbl))

ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                "ABcomp_PMW366-E1.E2_LFC_padj",
                                padjVal,"_lfc", localLFC,".pdf"),
                plot=p2, device="pdf",width=29,height=19,units="cm")





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
                            pattern="_merge_2000\\.E1\\.vecs\\.bw")
  E2files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.E2\\.vecs\\.bw")
  pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000.E1.vecs.bw")))]
  pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000.E2.vecs.bw")))]
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
  #### First eigen vector-----

  pcaSource="same"
  ######### 366 first Eigenvector ------
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


  ########## 366 second eigenvector ------
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


  #sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
  #                padj=padjVal,lfc=lfcVal,direction="both")
  #sigList<-lapply(listgr, as.data.frame)
  #sigList<-lapply(sigList, "[", ,c("compartment","log2FoldChange"))



  # upregulated
  sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                  padj=padjVal,lfc=lfcVal,direction="gt")

  sigList<-lapply(sigList, "[", ,c("compartment","log2FoldChange"))
  for(g in names(sigList)){ sigList[[g]]$SMC<-g }
  sigList<-do.call(rbind,sigList)
  sigList<-sigList[!is.na(sigList$compartment),]
  sigTbl<-sigList
  sigTbl$updown<-"up"


  # downregulated
  sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                  padj=padjVal, lfc= -lfcVal, direction="lt")

  sigList<-lapply(sigList, "[", ,c("compartment","log2FoldChange"))
  for(g in names(sigList)){ sigList[[g]]$SMC<-g }
  sigList<-do.call(rbind,sigList)
  sigList<-sigList[!is.na(sigList$compartment),]
  sigList$updown<-"down"
  sigTbl<-rbind(sigTbl,sigList)
  sigTbl$compartment<-as.factor(sigTbl$compartment)
  sigTbl$updown<-factor(sigTbl$updown, levels=c("up","down"))

  yminmax=c(0,max(abs(sigTbl$log2FoldChange)))
  p2<-ggplot(sigTbl,aes(x=compartment,y=abs(log2FoldChange),col=updown,fill=updown)) +
    geom_boxplot(notch=T, varwidth=T, position=position_dodge2(padding=0.2),outlier.size=0.4,
                 outlier.color="grey50") +
    facet_grid(cols=vars(SMC)) +
    ggtitle("Significantly changed genes by compartment") +
    theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
    theme(legend.title = element_text(size=10)) +
    scale_y_continuous(limits = yminmax) +
    scale_color_grey(start=0.2,end=0.2,guide="none")

  yminmax=c(0,median(abs(sigTbl$log2FoldChange))+quantile(abs(sigTbl$log2FoldChange))[4]*2)
  p3<-ggplot(sigTbl,aes(x=compartment,y=abs(log2FoldChange),col=updown,fill=updown)) +
    geom_boxplot(notch=T, varwidth=T, position=position_dodge2(padding=0.2),
                 outlier.shape=NA) +
    facet_grid(cols=vars(SMC)) +
    ggtitle("Significantly changed genes by compartment") +
    theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
    theme(legend.title = element_text(size=10)) +
    scale_y_continuous(limits = yminmax) +
    scale_color_grey(start=0.2,end=0.2,guide="none")


  p<-ggpubr::ggarrange(p2,p3,ncol=2,nrow=1)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                  "ABcomp_LFC_padj",
                                  padjVal,"_lfc", lfcVal,".pdf"),
                  plot=p, device="pdf",width=29,height=16,units="cm")



  # AB compartments - switching ---------------------------------------------

  ####
  ## sample specific compartments - changes between TEVonly and cs
  ####

  pcas<-data.frame(SMC=varOIlevels,
                   file=list.files(paste0(outPath,"/otherData"),
                                   pattern="_5000_laminDamID_pca2.bw"))
  listgr<-NULL
  for (grp in useContrasts){
    #grp=useContrasts[1]
    salmon<-readRDS(file=paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],
                                "_DESeq2_fullResults_p",padjVal,".rds"))
    pca2<-import.bw(paste0(outPath,"/otherData/",pcas$file[pcas$SMC==grp]))
    pca2control<-import.bw(paste0(outPath,"/otherData/",pcas$file[pcas$SMC==controlGrp]))

    salmon<-salmon[!is.na(salmon$chr),]
    salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

    salmongr<-sort(salmongr)


    salmongr<-assignGRtoAB(salmongr,pca2control,grName=controlGrp,pcaName=controlGrp)
    idx<-which(colnames(mcols(salmongr)) %in% c("pcaScore","compartment"))
    colnames(mcols(salmongr))[idx]<-paste(colnames(mcols(salmongr))[idx],"control",sep="_")
    salmongr<-assignGRtoAB(salmongr,pca2,grName=grp,pcaName=grp)
    salmongr$switch<-factor(paste0(salmongr$compartment_control,salmongr$compartment),levels=c("AA","BB","AB","BA"))
    listgr[[prettyGeneName(grp)]]<-salmongr
  }


  pairedCols<-c(brewer.pal(4,"Paired"))

  pdf(file=paste0(paste0(outPath,"/plots/",outputNamePrefix,
                         "ABcompSwitch_geneCount_padj",
                         padjVal,"_lfc", lfcVal,".pdf")),
      width=19, height=29, paper="a4")


  par(mfrow=c(3,1))
  # genes that change significantly
  sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                  padj=padjVal,lfc=lfcVal,direction="both")

  compartmentTable<-do.call(rbind,lapply(lapply(sigList, "[", ,"switch"),table))

  yminmax=c(0,max(compartmentTable))
  xx<-barplot(t(compartmentTable),beside=T,col=pairedCols,
              main="Significantly changed genes by compartment",cex.axis=1.2,
              cex.names=1.5, ylim=yminmax*1.1)
  legend("topright",legend = colnames(compartmentTable),fill=pairedCols)
  text(x=xx, y=t(compartmentTable), label=t(compartmentTable), pos=3,cex=1.1,col="black")

  par(mfrow=c(4,2))
  # upregulated genes
  sigListUp<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                    padj=padjVal,lfc=lfcVal,direction="gt")

  compartmentTableUp<-do.call(rbind,lapply(lapply(sigListUp, "[", ,"switch"),table))
  colnames(compartmentTableUp)<-paste0(colnames(compartmentTableUp),"_up")


  # downregulated genes
  sigListDown<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                      padj=padjVal, lfc= -lfcVal, direction="lt")

  compartmentTableDown<-do.call(rbind,lapply(lapply(sigListDown, "[", ,"switch"),table))
  colnames(compartmentTableDown)<-paste0(colnames(compartmentTableDown),"_down")

  compartmentTable<-cbind(compartmentTableUp,compartmentTableDown)

  yminmax=c(0,max(compartmentTable[,grep("AA|BB",colnames(compartmentTable))]))
  Acomp<-compartmentTable[,grep("AA",colnames(compartmentTable))]
  xx<-barplot(t(Acomp),beside=T,col=c("grey80","grey20"),
              main="Number of up/down regulated in AA compartment",cex.axis=1.2,
              cex.names=1.5, ylim=c(yminmax)*1.2)
  legend("top", legend=gsub("AA_","",colnames(Acomp)), fill=c("grey80","grey20"))
  text(x=xx, y=t(Acomp), label=t(Acomp), pos=3,cex=1.3,col="black")

  xx<-barplot(t(Acomp/rowSums(Acomp)),beside=F,col=c("grey80","grey20"),
              main="Fraction of up/down regulated in AA compartment",cex.axis=1.2,
              cex.names=1.5,space=0.8,ylim=c(0,1.1),bty='L')
  legend("bottomright", legend=gsub("AA_","",colnames(Acomp)), fill=c("grey80","grey20"),xpd=T)
  text(x=xx, y=0.97, label=t(rowSums(Acomp)), pos=3,cex=1.3,col="black")


  Bcomp<-compartmentTable[,grep("BB",colnames(compartmentTable))]
  xx<-barplot(t(Bcomp),beside=T,col=c("grey80","grey20"),
              main="Number of up/down regulated in BB compartment",cex.axis=1.2,
              cex.names=1.5, ylim=c(yminmax)*1.2)
  legend("top", legend=gsub("BB_","",colnames(Bcomp)), fill=c("grey80","grey20"))
  text(x=xx, y=t(Bcomp), label=t(Bcomp), pos=3,cex=1.3,col="black")

  xx<-barplot(t(Bcomp/rowSums(Bcomp)),beside=F,col=c("grey80","grey20"),
              main="Fraction of up/down regulated in BB compartment",cex.axis=1.2,
              cex.names=1.5,space=0.8,ylim=c(0,1.1),bty='L')
  legend("bottomright", legend=gsub("BB_","",colnames(Bcomp)), fill=c("grey80","grey20"),xpd=T)
  text(x=xx, y=0.97, label=t(rowSums(Bcomp)), pos=3,cex=1.3,col="black")


  yminmax=c(0,max(compartmentTable[,grep("AB|BA",colnames(compartmentTable))]))
  ABcomp<-compartmentTable[,grep("AB",colnames(compartmentTable))]
  xx<-barplot(t(ABcomp),beside=T,col=c("grey80","grey20"),
              main="Number of up/down regulated in AB compartment",cex.axis=1.2,
              cex.names=1.5, ylim=c(yminmax)*1.2)
  legend("top", legend=gsub("AB_","",colnames(ABcomp)), fill=c("grey80","grey20"))
  text(x=xx, y=t(ABcomp), label=t(ABcomp), pos=3,cex=1.3,col="black")

  xx<-barplot(t(ABcomp/rowSums(ABcomp)),beside=F,col=c("grey80","grey20"),
              main="Fraction of up/down regulated in AB compartment",cex.axis=1.2,
              cex.names=1.5,space=0.8,ylim=c(0,1.1),bty='L')
  legend("bottomright", legend=gsub("AB_","",colnames(ABcomp)), fill=c("grey80","grey20"),xpd=T)
  text(x=xx, y=0.97, label=t(rowSums(ABcomp)), pos=3,cex=1.3,col="black")


  BAcomp<-compartmentTable[,grep("BA",colnames(compartmentTable))]
  xx<-barplot(t(BAcomp),beside=T,col=c("grey80","grey20"),
              main="Number of up/down regulated in BA compartment",cex.axis=1.2,
              cex.names=1.5, ylim=c(yminmax)*1.2)
  legend("top", legend=gsub("BA_","",colnames(BAcomp)), fill=c("grey80","grey20"))
  text(x=xx, y=t(BAcomp), label=t(BAcomp), pos=3,cex=1.3,col="black")

  xx<-barplot(t(BAcomp/rowSums(BAcomp)),beside=F,col=c("grey80","grey20"),
              main="Fraction of up/down regulated in BA compartment",cex.axis=1.2,
              cex.names=1.5,space=0.8,ylim=c(0,1.1),bty='L')
  legend("bottomright", legend=gsub("BA_","",colnames(BAcomp)), fill=c("grey80","grey20"),xpd=T)
  text(x=xx, y=0.97, label=t(rowSums(BAcomp)), pos=3,cex=1.3,col="black")


  dev.off()






  ################-
  ## AB compartment by chromosome - switching between TEVonly and cs -----
  #################-

  # genes that change significantly
  sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                  padj=padjVal,lfc=lfcVal,direction="both")

  # count genes by category (chr & A/B)
  dfl<-lapply(sigList, function(x){x%>% dplyr::group_by(seqnames,switch,.drop=F) %>% tally()})

  # add name of SMC protein
  dfl<-do.call(rbind, mapply(cbind,dfl,"SMC"=names(dfl),SIMPLIFY=F))
  dfl$seqnames<-gsub("chr","",dfl$seqnames)
  ymax=max(dfl$n)
  p1<-ggplot(dfl,aes(x=seqnames,y=n,group=switch)) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=switch)) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_manual(values=pairedCols) +
    theme(legend.title = element_text(size=10)) +
    xlab("chr")+ylab("Number of genes") +
    ggtitle("Significantly changed genes per chromosome by compartment")


  dfl<-dfl[! (dfl$switch %in% c("AA","BB")),]
  dfl$switch<-droplevels(dfl$switch)
  ymax1=max(dfl$n)
  p1a<-ggplot(dfl,aes(x=seqnames,y=n,group=switch)) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=switch)) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_manual(values=pairedCols[3:4]) +
    theme(legend.title = element_text(size=10)) +
    xlab("chr")+ylab("Number of genes") +
    ggtitle("Significantly changed genes per chromosome by compartment")

  p<-ggpubr::ggarrange(p1,p1a,ncol=1,nrow=3)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                  "ABcompSwitch_countsPerChr_ABBA_padj",
                                  padjVal,"_lfc", lfcVal,".pdf"),
                  plot=p, device="pdf",width=19,height=29,units="cm")


  # upregulated genes
  sigListUp<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                    padj=padjVal,lfc=lfcVal,direction="gt")
  # count genes by category (chr & A/B)
  dflUp<-lapply(sigListUp, function(x){x%>% dplyr::group_by(seqnames,switch,.drop=F) %>% tally()})
  # add name of SMC protein
  dflUp<-do.call(rbind, mapply(cbind,dflUp,"SMC"=names(dflUp),SIMPLIFY=F))
  dflUp$seqnames<-gsub("chr","",dflUp$seqnames)
  dflUp$expression<-"up"

  # downregulated genes
  sigListDown<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                      padj=padjVal, lfc= -lfcVal, direction="lt")
  dflDown<-lapply(sigListDown, function(x){x%>% dplyr::group_by(seqnames,switch,.drop=F) %>% tally()})
  # add name of SMC protein
  dflDown<-do.call(rbind, mapply(cbind,dflDown,"SMC"=names(dflDown),SIMPLIFY=F))
  dflDown$seqnames<-gsub("chr","",dflDown$seqnames)
  dflDown$expression<-"down"

  dfl<-rbind(dflUp,dflDown)
  dfl$expression<-factor(dfl$expression,levels=c("up","down"))


  yminmax=c(0,max(dfl$n))
  p2<-ggplot(dfl[dfl$switch=="AA",],aes(x=seqnames,y=n,group=expression)) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=expression)) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
    theme(legend.title = element_text(size=10)) +
    xlab("chr")+ylab("Number of genes") + ylim(yminmax) +
    ggtitle("Up/down regulated genes in AA compartment")

  p3<-ggplot(dfl[dfl$switch=="BB",],aes(x=seqnames,y=n,group=expression)) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=expression)) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
    theme(legend.title = element_text(size=10)) +
    xlab("chr")+ylab("Number of genes") + ylim(yminmax) +
    ggtitle("Up/down regulated genes in BB compartment")

  yminmax=c(0,max(dfl$n[! (dfl$switch %in% c("AA","BB"))]))
  p4<-ggplot(dfl[dfl$switch=="AB",],aes(x=seqnames,y=n,group=expression)) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=expression)) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
    theme(legend.title = element_text(size=10)) +
    xlab("chr")+ylab("Number of genes") + ylim(yminmax) +
    ggtitle("Up/down regulated genes in AB compartment")

  p5<-ggplot(dfl[dfl$switch=="BA",],aes(x=seqnames,y=n,group=expression)) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=expression)) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
    theme(legend.title = element_text(size=10)) +
    xlab("chr")+ylab("Number of genes") + ylim(yminmax) +
    ggtitle("Up/down regulated genes in BA compartment")




  p<-ggpubr::ggarrange(p2,p3,p4,p5,ncol=2,nrow=2)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                  "ABcompSwitch_updownByChr_padj",
                                  padjVal,"_lfc", lfcVal,".pdf"),
                  plot=p, device="pdf",width=29,height=19,units="cm")

  ####
  ## AB comp LFC - switching between TEVonly and cs -----
  ####


  # upregulated
  sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                  padj=padjVal, lfc=lfcVal, direction="gt")

  sigList<-lapply(sigList, "[", ,c("switch","log2FoldChange"))
  for(g in names(sigList)){ sigList[[g]]$SMC<-g }
  sigList<-do.call(rbind,sigList)
  sigList<-sigList[!is.na(sigList$switch),]
  sigTbl<-sigList
  sigTbl$updown<-"up"


  # downregulated
  sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                  padj=padjVal, lfc= -lfcVal, direction="lt")

  sigList<-lapply(sigList, "[", ,c("switch","log2FoldChange"))
  for(g in names(sigList)){ sigList[[g]]$SMC<-g }
  sigList<-do.call(rbind,sigList)
  sigList<-sigList[!is.na(sigList$switch),]
  sigList$updown<-"down"
  sigTbl<-rbind(sigTbl,sigList)
  sigTbl$switch<-as.factor(sigTbl$switch)
  sigTbl$updown<-factor(sigTbl$updown, levels=c("up","down"))


  yminmax=c(0,max(abs(sigTbl$log2FoldChange)))
  p1<-ggplot(sigTbl,aes(x=switch,y=abs(log2FoldChange),col=updown,fill=updown)) +
    geom_boxplot(notch=T, varwidth=T, position=position_dodge2(padding=0.2),
                 outlier.size=0.4) +
    facet_grid(cols=vars(SMC)) +
    ggtitle("Log2 fold change by compartment") +
    theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
    theme(legend.title = element_text(size=10)) +
    scale_y_continuous(limits = yminmax) +
    scale_color_grey(start=0.7,end=0.3,guide="none")


  sigTbl<-sigTbl[! (sigTbl$switch %in% c("AA","BB")),]
  sigTbl$switch<-droplevels(sigTbl$switch)
  yminmax=c(0,max(abs(sigTbl$log2FoldChange)))
  p2<-ggplot(sigTbl,aes(x=switch,y=abs(log2FoldChange),col=updown,fill=updown)) +
    geom_boxplot(notch=T, varwidth=T, position=position_dodge2(padding=0.2),
                 outlier.shape=NA) +
    facet_grid(cols=vars(SMC)) +
    ggtitle("Log2 fold change by compartment") +
    theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
    theme(legend.title = element_text(size=10)) +
    scale_y_continuous(limits = yminmax) +
    scale_color_grey(start=0.7,end=0.3,guide="none")



  p<-ggpubr::ggarrange(p1,p2,ncol=2,nrow=1)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                  "ABcompSwitch_LFC_padj",
                                  padjVal,"_lfc", lfcVal,".pdf"),
                  plot=p, device="pdf",width=29,height=16,units="cm")


}

###########################-
# compartments - digitized ----------------------------------------------------
###########################-

RNAseqAndHiCsubset=c("aux_sdc3BG","dpy26","kle2","scc1","coh1")

if(all(RNAseqAndHiCsubset %in% useContrasts)){
  pcas<-data.frame(SMC=c("TEVonly","aux_sdc3BG","dpy26","kle2","scc1","coh1"),
                   strain =c("366","822","382","775","784","828"),
                   E1=NA, E2=NA)



  E1files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.saddle_cis_E1\\.digitized\\.tsv")
  E2files=list.files(paste0(outPath,"/otherData"),
                     pattern="_merge_2000\\.saddle_cis_E2\\.digitized\\.tsv")
  pcas$E1<-E1files[match(pcas$strain,unlist(strsplit(E1files,"_merge_2000\\.saddle_cis_E1\\.digitized\\.tsv")))]
  pcas$E2<-E2files[match(pcas$strain,unlist(strsplit(E2files,"_merge_2000\\.saddle_cis_E2\\.digitized\\.tsv")))]
  listgr<-NULL
  for (grp in RNAseqAndHiCsubset){
    #grp=useContrasts[1]
    salmon<-readRDS(file=paste0(outPath,"/rds/",fileNamePrefix,
                                contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
    pca1<-read.delim(paste0(outPath,"/otherData/",pcas$E1[pcas$SMC==grp]))
    pca2<-import.bw(paste0(outPath,"/otherData/",pcas$E2[pcas$SMC==grp]))

    salmon<-salmon[!is.na(salmon$chr),]
    salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

    salmongr<-sort(salmongr)

    salmongr<-assignGRtoAB(salmongr,pca1,grName=grp,pcaName=paste0("E1"))
    salmongr<-assignGRtoAB(salmongr,pca2,grName=grp,pcaName=paste0("E2"))
    listgr[[grp]]<-salmongr
  }







# seqplots heatmaps and averages ---------------------------------------------

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
#####################


######################-
# TADs ----
######################-
# loops<-rtracklayer::import(paste0(outPath,"/otherData/N2.allValidPairs.hic.5-10kbLoops.bedpe"), format="bedpe")
# head(loops)
# #extract the separate anchors
# grl<-zipup(loops)
# anchor1<-unlist(grl)[seq(1,2*length(grl),2)]
# anchor2<-unlist(grl)[seq(2,2*length(grl),2)]

# new SIP loops
loops<-read.delim(paste0(outPath,"/otherData/10kbLoops.txt"),header=T)
anchor1<-GRanges(seqnames=loops$chromosome1,ranges=IRanges(start=loops$x1,end=loops$x2-1),
                 strand="+")
anchor2<-GRanges(seqnames=loops$chromosome2,ranges=IRanges(start=loops$y1,end=loops$y2-1),
                 strand="-")
mcols(anchor1)<-loops[,c(8:15)]
mcols(anchor2)<-loops[,c(8:15)]
anchor1$loopNum<-paste0("loop",1:length(anchor1))
anchor2$loopNum<-paste0("loop",1:length(anchor2))

#make TADs
tads_in<-GRanges(seqnames=seqnames(anchor1),IRanges(start=end(anchor1)+1,end=start(anchor2)-1))

tads_in<-tads_in[width(tads_in)>100000]

#tads<-reduce(tads)
xtads<-tads_in[seqnames(tads_in)=="chrX"]
atads<-tads_in[seqnames(tads_in)!="chrX"]

flankSize<-200000

smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
                  useContrasts,"_lfc.bw")
names(smcRNAseq)<-useContrasts

pdf(file=paste0(outPath,"/plots/",outputNamePrefix,"tads-chrX_flank",
                  flankSize/1000,"kb.pdf"), width=11,
      height=9, paper="a4r")


p<-getPlotSetArray(tracks=c(smcRNAseq),
                   features=c(xtads),
                   refgenome="ce11", bin=10000L, xmin=flankSize,
                   xmax=flankSize, type="af",
                   xanchored=median(width(xtads)))



dd<-plotHeatmap(p,plotz=F)
heatmapQuantiles<-sapply(dd$HLST,quantile,c(0.05,0.95),na.rm=T)
roworder<-rev(order(lapply(dd$HLST,rowSums,na.rm=T)$X1))
minVal<-min(heatmapQuantiles[1,])
maxVal<-max(heatmapQuantiles[2,])
#layout(matrix(c(1), nrow = 2, ncol = 1, byrow = TRUE))
plotAverage(p,ylim=c(-0.5,1),main="ChrX TADs",
            error.estimates=ifelse(length(useContrasts>3),F,T))
plotHeatmap(p,main="ChrX TADs", plotScale="no", sortrows=T,
            clusters=1L,autoscale=F,zmin=minVal, zmax=maxVal,
            indi=F, sort_mids=T,sort_by=c(T,F,F),
            clspace=c("#00008B", "#FFFFE0","#8B0000"))


# reduced tads
p<-getPlotSetArray(tracks=c(smcRNAseq),
                   features=c(reduce(xtads)),
                   refgenome="ce11", bin=10000L, xmin=flankSize,
                   xmax=flankSize, type="af",
                   xanchored=median(width(xtads)))


dd<-plotHeatmap(p,plotz=F)
heatmapQuantiles<-sapply(dd$HLST,quantile,c(0.05,0.95),na.rm=T)
roworder<-rev(order(lapply(dd$HLST,rowSums,na.rm=T)$X1))
minVal<-min(heatmapQuantiles[1,])
maxVal<-max(heatmapQuantiles[2,])
#layout(matrix(c(1), nrow = 2, ncol = 1, byrow = TRUE))
plotAverage(p,ylim=c(-0.5,1),main="ChrX TADs (reduced)",error.estimates=F)
plotHeatmap(p,main="ChrX TADs (reduced)", plotScale="no", sortrows=T,
            clusters=1L,autoscale=F,zmin=minVal, zmax=maxVal,
            indi=F, sort_mids=T,sort_by=c(T,F,F),
            clspace=c("#00008B", "#FFFFE0","#8B0000"))


dev.off()




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


# ## old loops
# loops<-rtracklayer::import(paste0(outPath,"/otherData/N2.allValidPairs.hic.5-10kbLoops.bedpe"),
#                            format="bedpe")
# head(loops)
# #extract the separate anchors
# grl<-zipup(loops)
# anchor1<-unlist(grl)[seq(1,2*length(grl),2)]
# anchor2<-unlist(grl)[seq(2,2*length(grl),2)]


### new SIP loops
loops<-read.delim(paste0(outPath,"/otherData/10kbLoops.txt"),header=T)
anchor1<-GRanges(seqnames=loops$chromosome1,ranges=IRanges(start=loops$x1,end=loops$x2-1))
anchor2<-GRanges(seqnames=loops$chromosome2,ranges=IRanges(start=loops$y1,end=loops$y2-1))
mcols(anchor1)<-loops[,c(8:15)]
mcols(anchor2)<-loops[,c(8:15)]
anchor1$loopNum<-paste0("loop",1:length(anchor1))
anchor2$loopNum<-paste0("loop",1:length(anchor2))

# anchors<-c(anchor1,anchor2)
# #anchors<-reduce(anchors,min.gapwidth=0L)
# seqlevels(anchors)<-seqlevels(Celegans)[1:6]


#make TADs
tads<-GRanges(seqnames=seqnames(anchor1),IRanges(start=start(anchor1),end=end(anchor2)))
head(tads)
sort(width(tads))
tads<-reduce(tads)
sort(width(tads))

# find regions not in tads
notads<-gaps(tads)
sort(width(notads))


plotList<-list()
#grp=useContrasts[3]
for (grp in useContrasts){
  salmon<-readRDS(file=paste0(paste0(outPath,"/rds/",fileNamePrefix,
                                     contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))

  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

  salmongr<-sort(salmongr)

  ol<-findOverlaps(salmongr,tads,type="within")
  genesInTads<-salmongr[queryHits(ol)]

  ol<-findOverlaps(salmongr,notads)
  genesNotTads<-salmongr[queryHits(ol)]

  genesInTads$TADs<-"inside"
  genesNotTads$TADs<-"outside"
  df<-data.frame(c(genesInTads,genesNotTads))
  df<-df%>%dplyr::group_by(seqnames,TADs)%>%dplyr::mutate(count=n())

  plotList[[grp]]<-ggplot(df,aes(x=TADs,y=log2FoldChange,fill=TADs))+
    geom_boxplot(notch=T,outlier.shape=NA,varwidth=T)+
    facet_grid(.~seqnames) +ylim(c(-1,1))+
    ggtitle(grp)
}
p<-gridExtra::marrangeGrob(plotList,ncol=1,nrow=3)
ggsave(paste0(paste0(outPath,"/plots/",outputNamePrefix,"TADSinout_",
                     padjVal,".pdf")),
       width=9, height=11, paper="a4",plot=p,device="pdf")



#separate anchors from inside tads
tads_in<-reduce(GRanges(seqnames=seqnames(anchor1),IRanges(start=end(anchor1)+1,end=start(anchor2)-1)))
tads_in<-resize(tads_in,width=width(tads_in)-20000,fix="center")
anchors<-reduce(sort(c(anchor1,anchor2)))
#anchors<-resize(anchors,width=width(anchors)+20000,fix="center")
ol<-findOverlaps(anchors,tads_in)
anchors<-anchors[-queryHits(ol)]

plotList<-list()
#grp=useContrasts[3]
for (grp in useContrasts){
  salmon<-readRDS(file=paste0(paste0(outPath,"/rds/",fileNamePrefix,
                                     contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))

  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

  salmongr<-sort(salmongr)

  ol<-findOverlaps(salmongr,tads_in,type="within")
  insideTads<-salmongr[queryHits(ol)]

  ol<-findOverlaps(salmongr,anchors)
  atAnchors<-salmongr[queryHits(ol)]

  insideTads$TADs<-"TAD"
  atAnchors$TADs<-"Anchor"
  df<-data.frame(c(insideTads,atAnchors))
  df<-df%>%dplyr::group_by(seqnames,TADs)%>%dplyr::mutate(count=n())

  plotList[[grp]]<-ggplot(df,aes(x=TADs,y=log2FoldChange,fill=TADs))+
    geom_boxplot(notch=T,outlier.shape=NA,varwidth=T)+
    facet_grid(.~seqnames) +ylim(c(-1,1))+
    ggtitle(grp)
}
p<-gridExtra::marrangeGrob(plotList,ncol=1,nrow=3)
ggsave(paste0(paste0(outPath,"/plots/",outputNamePrefix,"TADSvAnchors_",
                     padjVal,".pdf")),
       width=9, height=11, paper="a4",plot=p,device="pdf")



####################-
## new SIP loops -----
####################-

# loops<-read.delim(paste0(outPath,"/otherData/10kbLoops.txt"),header=T)
# gr1<-GRanges(seqnames=loops$chromosome1,ranges=IRanges(start=loops$x1,end=loops$x2-1))
# gr2<-GRanges(seqnames=loops$chromosome2,ranges=IRanges(start=loops$y1,end=loops$y2-1))
# mcols(gr1)<-loops[,c(8:15)]
# mcols(gr2)<-loops[,c(8:15)]
# gr1$loopNum<-paste0("loop",1:length(gr1))
# gr2$loopNum<-paste0("loop",1:length(gr2))
#
# anchors<-c(gr1,gr2)
# #anchors<-reduce(anchors,min.gapwidth=0L)
# seqlevels(anchors)<-seqlevels(Celegans)[1:6]
#
# smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
#                   useContrasts,"_lfc.bw")
# names(smcRNAseq)<-useContrasts
# for(grp in useContrasts){
#   rnaSeq<-import.bw(smcRNAseq[[grp]])
#   cov<-coverage(rnaSeq,weight="score")
#   anchors<-binnedAverage(anchors,cov,grp)
# }
#
#
# pdf(paste0(outPath, "/plots/",outputNamePrefix,"RNAseq_loopsMetrics.pdf"),
#            width=11,height=8,paper="a4r")
# winSize=10000
# metricsOI<-c("APScoreAvg","ProbabilityofEnrichment","RegAPScoreAvg","Avg_diffMaxNeihgboor_1", "Avg_diffMaxNeihgboor_2","avg","std","value")
# plotList<-list()
# for (colOI in metricsOI){
#   df<-data.frame(anchors)
#   colnames(df)<-c(colnames(df)[1:5],colnames(mcols(anchors)))
#   df<-tidyr::pivot_longer(df,useContrasts,names_to="SMC",values_to="LFC")
#   plotList[[colOI]]<-ggplot(df,aes_string(x=colOI,y="LFC",col="SMC")) + geom_point() +
#     ylab("Average RNAseq score (10kb window)") + xlab(colOI) +
#     ggtitle(paste0("Average RNAseq in ",winSize/1000,"kb window vs ",colOI))
# }
# p<-ggpubr::ggarrange(plotlist=plotList,ncol=2,nrow=2)
# print(p)
# dev.off()



###################-
## binned signal around anchors
###################-

# new SIP loops
loops<-read.delim(paste0(outPath,"/otherData/10kbLoops.txt"),header=T)
anchor1<-GRanges(seqnames=loops$chromosome1,ranges=IRanges(start=loops$x1,end=loops$x2-1),
             strand="+")
anchor2<-GRanges(seqnames=loops$chromosome2,ranges=IRanges(start=loops$y1,end=loops$y2-1),
             strand="-")
mcols(anchor1)<-loops[,c(8:15)]
mcols(anchor2)<-loops[,c(8:15)]
anchor1$loopNum<-paste0("loop",1:length(anchor1))
anchor2$loopNum<-paste0("loop",1:length(anchor2))

anchors<-c(anchor1,anchor2)
anchors<-reduce(anchors,min.gapwidth=0L)
seqlevels(anchors)<-seqlevels(Celegans)[1:6]

xanchors<-anchors[seqnames(anchors)=="chrX"]

smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
                  useContrasts,"_lfc.bw")
names(smcRNAseq)<-useContrasts

p<-avrSignalBins(xanchors, bwFiles=smcRNAseq, winSize=10000,numWins=10)
p<-p+ggtitle("Average RNAseq LFC around chrX SIP loop anchors")
print(p)
ggsave(paste0(outPath, "/plots/",outputNamePrefix,"binnedRNAseq_nearSIPanchors.pdf"), p,
       device="pdf",width=19,height=29,units="cm")



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




##########################-
########### moustache loops------
##########################-

### new Moustache loops
#mustacheBatch="PMW366"
mustacheBatch="PMW382"

loops<-import(paste0(outPath,"/otherData/",mustacheBatch,"_2k_mustache_filtered.bedpe"),format="bedpe")
grl<-zipup(loops)
anchor1<-do.call(c,lapply(grl,"[",1))
anchor2<-do.call(c,lapply(grl,"[",2))
mcols(anchor1)<-mcols(loops)
mcols(anchor2)<-mcols(loops)

anchor1$loopNum<-paste0("loop",1:length(anchor1))
anchor2$loopNum<-paste0("loop",1:length(anchor2))

# anchors<-c(anchor1,anchor2)
# #anchors<-reduce(anchors,min.gapwidth=0L)
# seqlevels(anchors)<-seqlevels(Celegans)[1:6]

#make TADs
tads<-GRanges(seqnames=seqnames(anchor1),IRanges(start=start(anchor1),end=end(anchor2)))
head(tads)
sort(width(tads))
tads<-reduce(tads)
sort(width(tads))

# find regions not in tads
notads<-gaps(tads)
sort(width(notads))


plotList<-list()
#grp=useContrasts[3]
for (grp in useContrasts){
  salmon<-readRDS(file=paste0(paste0(outPath,"/rds/",fileNamePrefix,
                                     contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))

  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

  salmongr<-sort(salmongr)

  ol<-findOverlaps(salmongr,tads,type="within")
  genesInTads<-salmongr[queryHits(ol)]

  ol<-findOverlaps(salmongr,notads)
  genesNotTads<-salmongr[queryHits(ol)]

  genesInTads$Loops<-"inside"
  genesNotTads$Loops<-"outside"
  df<-data.frame(c(genesInTads,genesNotTads))
  df<-df%>%dplyr::group_by(seqnames,Loops)%>%dplyr::mutate(count=n())

  plotList[[grp]]<-ggplot(df,aes(x=Loops,y=log2FoldChange,fill=Loops))+
    geom_boxplot(notch=T,outlier.shape=NA,varwidth=T)+
    facet_grid(.~seqnames) +ylim(c(-1,1))+
    ggtitle(grp)
}
p<-gridExtra::marrangeGrob(plotList,ncol=1,nrow=3)
ggsave(paste0(paste0(outPath,"/plots/",outputNamePrefix,"LoopsInOut_",mustacheBatch,"Mustache_",
                     padjVal,".pdf")),
       width=9, height=11, paper="a4",plot=p,device="pdf")



#separate anchors from inside tads
tads_in<-GRanges(seqnames=seqnames(anchor1),IRanges(start=end(anchor1)+1,end=start(anchor2)-1))
# tads_in<-reduce(tads_in)
tads_in<-resize(tads_in,width=width(tads_in)-20000,fix="center")
anchors<-sort(c(anchor1,anchor2))
anchors<-resize(anchors,width=20000,fix="center")
#anchors<-reduce(anchors)
#ol<-findOverlaps(anchors,tads_in)
#anchors<-anchors[-queryHits(ol)]

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

  ol<-findOverlaps(salmongr,tads_in,type="within")
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

ggsave(paste0(paste0(outPath,"/plots/",outputNamePrefix,"LoopsvAnchors_",mustacheBatch,"-Mostache_",
                     padjVal,".pdf")),
       width=9, height=11, paper="a4",plot=p,device="pdf")

## focus on chrX loops
dataTbl<-do.call(rbind,dataList)
xchr<-dataTbl[dataTbl$seqnames=="chrX",]

xchr$SMC<-factor(xchr$SMC,levels=useContrasts)
p1<-ggplot(xchr,aes(x=Loops,y=log2FoldChange,fill=Loops))+
  geom_boxplot(notch=T,outlier.shape=NA,varwidth=T)+
  facet_grid(~SMC) +ylim(c(-1,1))+
  ggtitle(paste0("LFC at ",mustacheBatch," anchors Vs inside loops in chrX")) +
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
ggsave(paste0(paste0(outPath,"/plots/",outputNamePrefix,"LoopsvAnchorsXchr_",mustacheBatch,"-Mostache_",
                     padjVal,".pdf")),
       width=11, height=5, paper="a4r",plot=p,device="pdf")



################-
## compare anchors
################-

loops366<-import(paste0(outPath,"/otherData/PMW366_2k_mustache_filtered.bedpe"),format="bedpe")

loops382<-import(paste0(outPath,"/otherData/PMW382_2k_mustache_filtered.bedpe"),format="bedpe")



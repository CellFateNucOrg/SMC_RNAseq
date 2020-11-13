library(rtracklayer)
#library(ggVennDiagram)
library(ggplot2)
#library(EnhancedVolcano)
library(zoo)
library(dplyr)
library(ggpubr)
library(genomation)

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
## AB compartments
############################

pca2<-import.bw("./otherData/N2_5000_DamID_pca2.bw")

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

lfcVal=0.5
padjVal=0.05
listgr<-NULL
for (grp in groupsOI){
  #grp=groupsOI[1]
  salmon<-readRDS(file=paste0("./rds/salmon_",grp,"_DESeq2_fullResults.rds"))

  sigGenes<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="both",
                        chr="all", nameChrCol="chr"))

  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

  salmongr<-sort(salmongr)

  ol<-as.data.frame(findOverlaps(salmongr,pca2,ignore.strand=T))
  ol$subjectScore<-pca2$score[ol$subjectHits]
  pcaScore<-ol %>% group_by(queryHits)%>% summarise(pcaScore=mean(subjectScore))

  #salmongr$pcaScore<-NA
  salmongr$pcaScore[pcaScore$queryHits]<-pcaScore$pcaScore
  salmongr$compartment<-as.factor(ifelse(salmongr$pcaScore>0,"B","A"))
  forBG<-salmongr[! is.na (salmongr$compartment)]
  forBG$score<-ifelse(forBG$compartment=="A",1,-1)
  export(forBG,con=paste0("./tracks/checkCompartments_",grp,".bedGraph"),
         format="bedGraph")
  listgr[[prettyGeneName(grp)]]<-salmongr
}



pdf(file=paste0(paste0(outPath,"/plots/ABcomp_geneCount_padj",
                       padjVal,"_lfc", lfcVal,".pdf")),
    width=19, height=4, paper="a4r")
par(mfrow=c(1,3))
# genes that change significantly
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="both")

compartmentTable<-do.call(rbind,lapply(lapply(sigList, "[", ,"compartment"),table))

yminmax=c(0,max(compartmentTable))
barplot(t(compartmentTable),beside=T,col=c("grey80","grey20"),
        main="Number of genes that change significantly",cex.axis=1.2,
        cex.names=1.5, ylim=yminmax)
legend("topleft",legend = colnames(compartmentTable),fill=c("grey80","grey20"))

# upregulated genes
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="gt")

compartmentTable<-do.call(rbind,lapply(lapply(sigList, "[", ,"compartment"),table))

barplot(t(compartmentTable),beside=T,col=c("grey80","grey20"),
        main="Upregulated genes by compartment",cex.axis=1.2,
        cex.names=1.5, ylim=yminmax)
legend("topleft", legend=colnames(compartmentTable), fill=c("grey80","grey20"))


# downregulated genes
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=-lfcVal,direction="lt")

compartmentTable<-do.call(rbind,lapply(lapply(sigList, "[", ,"compartment"),table))

barplot(t(compartmentTable),beside=T,col=c("grey80","grey20"),
        main="Downregulated genes by compartment",cex.axis=1.2,
        cex.names=1.5, ylim=yminmax)
legend("topleft", legend=colnames(compartmentTable), fill=c("grey80","grey20"))
dev.off()


# genes that change significantly
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="both")

sigList<-lapply(sigList, "[", ,c("compartment","log2FoldChange"))
#sigList$SMC<-NA
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList<-sigList[!is.na(sigList$compartment),]
sigList$compartment<-as.factor(sigList$compartment)

yminmax=max(abs(min(sigList$log2FoldChange)),max(sigList$log2FoldChange))
yminmax<-c(-yminmax,yminmax)
p1<-ggplot(sigList,aes(x=compartment,y=log2FoldChange,fill=compartment)) +
  geom_violin() + facet_grid(cols=vars(SMC)) +
  ylim(yminmax) +
  ggtitle("Significantly changed genes by compartment") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3)


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
                padj=padjVal,lfc=-lfcVal,direction="lt")

sigList<-lapply(sigList, "[", ,c("compartment","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList<-sigList[!is.na(sigList$compartment),]
sigList$updown<-"down"
sigTbl<-rbind(sigTbl,sigList)
sigTbl$compartment<-as.factor(sigTbl$compartment)
sigTbl$updown<-as.factor(sigTbl$updown)

yminmax=c(-max(abs(sigTbl$log2FoldChange)),max(abs(sigTbl$log2FoldChange)))
p2<-ggplot(sigTbl,aes(x=compartment,y=log2FoldChange,col=updown,fill=compartment)) +
  geom_boxplot(notch=T, varwidth=T, position=position_dodge(0),outlier.size=0.3) +
  facet_grid(cols=vars(SMC)) + ylim(yminmax) +
  ggtitle("Significantly changed genes by compartment") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide=F)

p<-ggpubr::ggarrange(p1,p2,ncol=2,nrow=1)
ggplot2::ggsave(filename=paste0(outPath, "/plots/ABcomp_LFC_padj",
                          padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=29,height=16,units="cm")



#####################
## anchors
#####################

loops<-import("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/2020_RNA_seq/otherData/N2.allValidPairs.hic.5-10kbLoops.bedpe",format="bedpe")


anchors<-unlist(zipup(loops))
anchors<-anchors[seqnames(anchors)=="chrX"]

fullSet<-reduce(anchors)

tenkb<-fullSet[width(fullSet)==10000]
fivekb<-fullSet[width(fullSet)==5000]



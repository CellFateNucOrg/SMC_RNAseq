library(rtracklayer)
#library(ggVennDiagram)
library(ggplot2)
#library(EnhancedVolcano)
library(zoo)
library(dplyr)
library(ggpubr)
library(genomation)
library(seqplots)
library(RColorBrewer)
source("functions.R")

outPath="."
padjVal=0.05
lfcVal=0.5
plotPDFs=F

fileList<-read.table(paste0(outPath,"/fastqList.txt"),stringsAsFactors=F,header=T)

# extract the strain variable
strain<-factor(as.character(unique(fileList$sampleName),levels=c("366","382","775","784")))
SMC<-strain
levels(SMC)<-c("TEVonly","dpy26cs","kle2cs","scc1cs")

controlGrp<-levels(SMC)[1] # control group
groupsOI<-levels(SMC)[-1]


# AB compartments - N2 ----------------------------------------------------


####
## AB compartments
####

####
## N2 compartments
####
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

listgr<-NULL
for (grp in groupsOI){
  #grp=groupsOI[1]
  salmon<-readRDS(file=paste0("./rds/salmon_",grp,"_DESeq2_fullResults.rds"))

  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

  salmongr<-sort(salmongr)

  salmongr<-assignGRtoAB(salmongr,pca2,grName=grp,pcaName="N2")
  listgr[[prettyGeneName(grp)]]<-salmongr
}



pdf(file=paste0(paste0(outPath,"/plots/ABcomp_N2_geneCount_padj",
                       padjVal,"_lfc", lfcVal,".pdf")),
    width=19, height=4, paper="a4r")
par(mfrow=c(1,3))
# genes that change significantly
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="both")

compartmentTable<-do.call(rbind,lapply(lapply(sigList, "[", ,"compartment"),table))

yminmax=c(0,max(compartmentTable))
barplot(t(compartmentTable),beside=T,col=c("grey80","grey20"),
        main="Significantly changed genes by N2 compartment",cex.axis=1.2,
        cex.names=1.5, ylim=yminmax)
legend("topleft",legend = colnames(compartmentTable),fill=c("grey80","grey20"))

# upregulated genes
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="gt")

compartmentTable<-do.call(rbind,lapply(lapply(sigList, "[", ,"compartment"),table))

barplot(t(compartmentTable),beside=T,col=c("grey80","grey20"),
        main="Upregulated genes by N2 compartment",cex.axis=1.2,
        cex.names=1.5, ylim=yminmax)
legend("topleft", legend=colnames(compartmentTable), fill=c("grey80","grey20"))


# downregulated genes
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=-lfcVal,direction="lt")

compartmentTable<-do.call(rbind,lapply(lapply(sigList, "[", ,"compartment"),table))

barplot(t(compartmentTable),beside=T,col=c("grey80","grey20"),
        main="Downregulated genes by N2 compartment",cex.axis=1.2,
        cex.names=1.5, ylim=yminmax)
legend("topleft", legend=colnames(compartmentTable), fill=c("grey80","grey20"))
dev.off()


####
## AB compartment by chromosome
####
# genes that change significantly
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="both")

# count genes by category (chr & A/B)
dfl<-lapply(sigList, function(x){x%>% group_by(seqnames,compartment) %>% tally()})
# add name of SMC protein
dfl<-do.call(rbind, mapply(cbind,dfl,"SMC"=names(dfl),SIMPLIFY=F))
dfl$seqnames<-gsub("chr","",dfl$seqnames)
ymax=max(dfl$n)
p1<-ggplot(dfl,aes(x=seqnames,y=n,group=compartment)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=compartment)) +
  facet_grid(cols=vars(SMC)) +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Significantly changed genes per chromosome by N2 compartment")


# upregulated genes
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="gt")


# count genes by category (chr & A/B)
dfl<-lapply(sigList, function(x){x%>% group_by(seqnames,compartment) %>% tally()})
# add name of SMC protein
dfl<-do.call(rbind, mapply(cbind,dfl,"SMC"=names(dfl),SIMPLIFY=F))
dfl$seqnames<-gsub("chr","",dfl$seqnames)
yminmax=c(0,max(dfl$n))
p2<-ggplot(dfl,aes(x=seqnames,y=n,group=compartment)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=compartment)) +
  facet_grid(cols=vars(SMC)) +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Upregulated genes per chromosome by N2 compartment")


# downregulated genes
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=-lfcVal,direction="lt")


dfl<-lapply(sigList, function(x){x%>% group_by(seqnames,compartment) %>% tally()})
# add name of SMC protein
dfl<-do.call(rbind, mapply(cbind,dfl,"SMC"=names(dfl),SIMPLIFY=F))
dfl$seqnames<-gsub("chr","",dfl$seqnames)
yminmax=c(0,max(dfl$n))
p3<-ggplot(dfl,aes(x=seqnames,y=n,group=compartment)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=compartment)) +
  facet_grid(cols=vars(SMC),switch="x") +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Downregulated genes per chromosome by N2 compartment") +
  scale_y_reverse(limits=c(ymax,0)) + scale_x_discrete(position = "top")

p<-ggpubr::ggarrange(p1,p2,p3,ncol=1,nrow=3)
ggplot2::ggsave(filename=paste0(outPath, "/plots/ABcomp_N2_countsPerChr_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=19,height=29,units="cm")



####
## AB comp LFC
####

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
  ggtitle("Significantly changed genes by N2 compartment") +
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
  ggtitle("Significantly changed genes by N2 compartment") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide=F)

p<-ggpubr::ggarrange(p1,p2,ncol=2,nrow=1)
ggplot2::ggsave(filename=paste0(outPath, "/plots/ABcomp_N2_LFC_padj",
                          padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=29,height=16,units="cm")






# AB compartments - sample specific ---------------------------------------

####
## sample specific compartments
####

pcas<-data.frame(SMC=SMC,
                 file=list.files("./otherData",
                                pattern="Illumina_5000.cool_DamID.pca2.bw"))

listgr<-NULL
for (grp in groupsOI){
  #grp=groupsOI[1]
  salmon<-readRDS(file=paste0("./rds/salmon_",grp,"_DESeq2_fullResults.rds"))
  pca2<-import.bw(paste0(outPath,"/otherData/",pcas$file[pcas$SMC==grp]))

  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

  salmongr<-sort(salmongr)

  salmongr<-assignGRtoAB(salmongr,pca2,grName=grp,pcaName=grp)
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
        main="Significantly changed genes by compartment",cex.axis=1.2,
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


####
## AB compartment by chromosome
####

# genes that change significantly
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="both")

# count genes by category (chr & A/B)
dfl<-lapply(sigList, function(x){x%>% group_by(seqnames,compartment) %>% tally()})
# add name of SMC protein
dfl<-do.call(rbind, mapply(cbind,dfl,"SMC"=names(dfl),SIMPLIFY=F))
dfl$seqnames<-gsub("chr","",dfl$seqnames)
ymax=max(dfl$n)
p1<-ggplot(dfl,aes(x=seqnames,y=n,group=compartment)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=compartment)) +
  facet_grid(cols=vars(SMC)) +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Significantly changed genes per chromosome by compartment")


# upregulated genes
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="gt")


# count genes by category (chr & A/B)
dfl<-lapply(sigList, function(x){x%>% group_by(seqnames,compartment) %>% tally()})
# add name of SMC protein
dfl<-do.call(rbind, mapply(cbind,dfl,"SMC"=names(dfl),SIMPLIFY=F))
dfl$seqnames<-gsub("chr","",dfl$seqnames)
yminmax=c(0,max(dfl$n))
p2<-ggplot(dfl,aes(x=seqnames,y=n,group=compartment)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=compartment)) +
  facet_grid(cols=vars(SMC)) +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Upregulated genes per chromosome by compartment")


# downregulated genes
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=-lfcVal,direction="lt")


dfl<-lapply(sigList, function(x){x%>% group_by(seqnames,compartment) %>% tally()})
# add name of SMC protein
dfl<-do.call(rbind, mapply(cbind,dfl,"SMC"=names(dfl),SIMPLIFY=F))
dfl$seqnames<-gsub("chr","",dfl$seqnames)
yminmax=c(0,max(dfl$n))
p3<-ggplot(dfl,aes(x=seqnames,y=n,group=compartment)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=compartment)) +
  facet_grid(cols=vars(SMC),switch="x") +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Downregulated genes per chromosome by compartment") +
  scale_y_reverse(limits=c(ymax,0)) + scale_x_discrete(position = "top")

p<-ggpubr::ggarrange(p1,p2,p3,ncol=1,nrow=3)
ggplot2::ggsave(filename=paste0(outPath, "/plots/ABcomp_countsPerChr_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=19,height=29,units="cm")



####
## AB comp LFC
####

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






# AB compartments - switching ---------------------------------------------

####
## sample specific compartments - changes between TEVonly and cs
####

pcas<-data.frame(SMC=SMC,
                 file=list.files("./otherData",
                                 pattern="Illumina_5000.cool_DamID.pca2.bw"))

listgr<-NULL
for (grp in groupsOI){
  #grp=groupsOI[1]
  salmon<-readRDS(file=paste0("./rds/salmon_",grp,"_DESeq2_fullResults.rds"))
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

pdf(file=paste0(paste0(outPath,"/plots/ABcompSwitch_geneCount_padj",
                       padjVal,"_lfc", lfcVal,".pdf")),
    width=19, height=4, paper="a4r")
par(mfrow=c(1,3))
# genes that change significantly
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="both")

compartmentTable<-do.call(rbind,lapply(lapply(sigList, "[", ,"switch"),table))

yminmax=c(0,max(compartmentTable))
barplot(t(compartmentTable),beside=T,col=pairedCols,
        main="Significantly changed genes by compartment",cex.axis=1.2,
        cex.names=1.5, ylim=yminmax)
legend("topright",legend = colnames(compartmentTable),fill=pairedCols)

# upregulated genes
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="gt")

compartmentTable<-do.call(rbind,lapply(lapply(sigList, "[", ,"switch"),table))

barplot(t(compartmentTable),beside=T,col=pairedCols,
        main="Upregulated genes by compartment",cex.axis=1.2,
        cex.names=1.5, ylim=yminmax)
legend("topright", legend=colnames(compartmentTable), fill=pairedCols)


# downregulated genes
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=-lfcVal,direction="lt")

compartmentTable<-do.call(rbind,lapply(lapply(sigList, "[", ,"switch"),table))

barplot(t(compartmentTable),beside=T,col=pairedCols,
        main="Downregulated genes by compartment",cex.axis=1.2,
        cex.names=1.5, ylim=yminmax)
legend("topright", legend=colnames(compartmentTable), fill=pairedCols)
dev.off()


################-
## AB compartment by chromosome - changes between TEVonly and cs
#################-

# genes that change significantly
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="both")

# count genes by category (chr & A/B)
dfl<-lapply(sigList, function(x){x%>% group_by(seqnames,switch,.drop=F) %>% tally()})

# add name of SMC protein
dfl<-do.call(rbind, mapply(cbind,dfl,"SMC"=names(dfl),SIMPLIFY=F))
dfl$seqnames<-gsub("chr","",dfl$seqnames)
ymax=max(dfl$n)
p1<-ggplot(dfl,aes(x=seqnames,y=n,group=switch)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=switch)) +
  facet_grid(cols=vars(SMC)) +
  theme_minimal() + scale_fill_manual(values=pairedCols) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Significantly changed genes per chromosome by compartment")

dfl<-dfl[! (dfl$switch %in% c("AA","BB")),]
dfl$switch<-droplevels(dfl$switch)
ymax1=max(dfl$n)
p1a<-ggplot(dfl,aes(x=seqnames,y=n,group=switch)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=switch)) +
  facet_grid(cols=vars(SMC)) +
  theme_minimal() + scale_fill_manual(values=pairedCols[3:4]) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Significantly changed genes per chromosome by compartment")

# upregulated genes
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="gt")


# count genes by category (chr & A/B)
dfl<-lapply(sigList, function(x){x%>% group_by(seqnames,switch,.drop=F) %>% tally()})
# add name of SMC protein
dfl<-do.call(rbind, mapply(cbind,dfl,"SMC"=names(dfl),SIMPLIFY=F))
dfl$seqnames<-gsub("chr","",dfl$seqnames)
yminmax=c(0,max(dfl$n))
p2<-ggplot(dfl,aes(x=seqnames,y=n,group=switch)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=switch)) +
  facet_grid(cols=vars(SMC)) +
  theme_minimal() + scale_fill_manual(values=pairedCols) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Upregulated genes per chromosome by compartment")

dfl<-dfl[! (dfl$switch %in% c("AA","BB")), ]
dfl$switch<-droplevels(dfl$switch)
ymax1=max(dfl$n)

p2a<-ggplot(dfl,aes(x=seqnames,y=n,group=switch)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=switch)) +
  facet_grid(cols=vars(SMC)) +
  theme_minimal() + scale_fill_manual(values=pairedCols[3:4]) +
  xlab("chr")+ylab("Number of genes")
  ggtitle("Upregulated genes per chromosome by compartment")

# downregulated genes
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=-lfcVal,direction="lt")


dfl<-lapply(sigList, function(x){x%>% group_by(seqnames,switch,.drop=F) %>% tally()})
# add name of SMC protein
dfl<-do.call(rbind, mapply(cbind,dfl,"SMC"=names(dfl),SIMPLIFY=F))
dfl$seqnames<-gsub("chr","",dfl$seqnames)
yminmax=c(0,max(dfl$n))
p3<-ggplot(dfl,aes(x=seqnames,y=n,group=switch)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=switch)) +
  facet_grid(cols=vars(SMC),switch="x") +
  theme_minimal() + scale_fill_manual(values=pairedCols) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Downregulated genes per chromosome by compartment") +
  scale_y_reverse(limits=c(ymax,0)) + scale_x_discrete(position = "top")

dfl<-dfl[! (dfl$switch %in% c("AA","BB")), ]
dfl$switch<-droplevels(dfl$switch)

p3a<-ggplot(dfl,aes(x=seqnames,y=n,group=switch)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=switch)) +
  facet_grid(cols=vars(SMC),switch="x") +
  theme_minimal() + scale_fill_manual(values=pairedCols[3:4]) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Downregulated genes per chromosome by compartment") +
  scale_y_reverse(limits=c(ymax1,0)) + scale_x_discrete(position = "top")


p<-ggpubr::ggarrange(p1,p2,p3,ncol=1,nrow=3)
ggplot2::ggsave(filename=paste0(outPath, "/plots/ABcompSwitch_countsPerChr_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=19,height=29,units="cm")


p<-ggpubr::ggarrange(p1a,p2a,p3a,ncol=1,nrow=3)
ggplot2::ggsave(filename=paste0(outPath, "/plots/ABcompSwitch_countsPerChr_ABBA_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=19,height=29,units="cm")


####
## AB comp LFC - changes between TEVonly and cs
####

# genes that change significantly
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="both")

sigList<-lapply(sigList, "[", ,c("switch","log2FoldChange"))
#sigList$SMC<-NA
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList<-sigList[!is.na(sigList$switch),]

yminmax=max(abs(min(sigList$log2FoldChange)),max(sigList$log2FoldChange))
yminmax<-c(-yminmax,yminmax)
p1<-ggplot(sigList,aes(x=switch,y=log2FoldChange,fill=switch)) +
  geom_violin() + facet_grid(cols=vars(SMC)) +
  ylim(yminmax) +
  ggtitle("Significantly changed genes by compartment") +
  theme_minimal() + scale_fill_manual(values=pairedCols)

sigList<-sigList[! (sigList$switch %in% c("AA","BB")),]
sigList$switch<-droplevels(sigList$switch)
yminmax1=max(abs(min(sigList$log2FoldChange)),max(sigList$log2FoldChange))
yminmax1<-c(-yminmax1,yminmax1)
p1a<-ggplot(sigList,aes(x=switch,y=log2FoldChange,fill=switch)) +
  geom_violin() + facet_grid(cols=vars(SMC)) +
  ylim(yminmax) +
  ggtitle("Significantly changed genes by compartment") +
  theme_minimal() + scale_fill_manual(values=pairedCols[3:4])


# upregulated
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="gt")

sigList<-lapply(sigList, "[", ,c("switch","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList<-sigList[!is.na(sigList$switch),]
sigTbl<-sigList
sigTbl$updown<-"up"


# downregulated
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=-lfcVal,direction="lt")

sigList<-lapply(sigList, "[", ,c("switch","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList<-sigList[!is.na(sigList$switch),]
sigList$updown<-"down"
sigTbl<-rbind(sigTbl,sigList)
#sigTbl$compartment<-as.factor(sigTbl$compartment)
sigTbl$updown<-as.factor(sigTbl$updown)

yminmax=c(-max(abs(sigTbl$log2FoldChange)),max(abs(sigTbl$log2FoldChange)))
p2<-ggplot(sigTbl,aes(x=switch,y=log2FoldChange,col=updown,fill=switch)) +
  geom_boxplot(notch=T, varwidth=T, position=position_dodge(0),outlier.size=0.3) +
  facet_grid(cols=vars(SMC)) + #ylim(yminmax) +
  ggtitle("Significantly changed genes by compartment") +
  theme_minimal() + scale_fill_manual(values=pairedCols) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide=F)

sigTbl<-sigTbl[! (sigTbl$switch %in% c("AA","BB")),]
sigTbl$switch<-droplevels(sigTbl$switch)
yminmax1=max(abs(min(sigTbl$log2FoldChange)),max(sigTbl$log2FoldChange))
yminmax1<-c(-yminmax1,yminmax1)

p2a<-ggplot(sigTbl,aes(x=switch,y=log2FoldChange,col=updown,fill=switch)) +
  geom_boxplot(notch=T, varwidth=T, position=position_dodge(0),outlier.size=0.3) +
  facet_grid(cols=vars(SMC)) + #ylim(yminmax1) +
  ggtitle("Significantly changed genes by compartment") +
  theme_minimal() + scale_fill_manual(values=pairedCols[3:4]) +
  scale_y_continuous(limits = yminmax1) +
  scale_color_grey(start=0.2,end=0.2,guide=F)



p<-ggpubr::ggarrange(p1,p2,ncol=2,nrow=1)
ggplot2::ggsave(filename=paste0(outPath, "/plots/ABcompSwitch_LFC_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=29,height=16,units="cm")

p<-ggpubr::ggarrange(p1a,p2a,ncol=2,nrow=1)
ggplot2::ggsave(filename=paste0(outPath, "/plots/ABcompSwitch_LFC_ABBA_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=29,height=16,units="cm")




# anchors - genomation ----------------------------------------------------

# ###
# ## anchors
# ###
#
# loops<-import(paste0(outPath,"/otherData/N2.allValidPairs.hic.5-10kbLoops.bedpe"),format="bedpe")
#
# anchors<-zipup(loops)
# anchors<-unlist(anchors)
# strand(anchors)[seq(1,length(anchors),by=2)]<-"+"
# strand(anchors)[seq(2,length(anchors),by=2)]<-"-"
#
# anchors<-reduce(sort(anchors))
# #anchors$region<-1:length(anchors)
# anchors$chr<-seqnames(anchors)
# keepAnchors<-anchors
# flankSize<-10000
# par(mfrow=c(3,1))
# for (grp in groupsOI){
#   anchors<-keepAnchors
#   smcRNAseq<-import(paste0(outPath,"/tracks/salmon_",grp,
#                            "_wt_lfc.bw"),
#                            format="bigwig")
#   pdf(file=paste0(outPath,"/plots/anchors_",grp,".pdf"),
#       width=11, height=29,paper="a4")
#   upstream<-trim(flank(anchors,width=flankSize))
#   downstream<-trim(flank(anchors,width=flankSize,start=F))
#
#   upstream<-upstream[order(upstream$region)]
#   downstream<-downstream[order(downstream$region)]
#
#   up2<-trim(flank(upstream,width=flankSize))
#   down2<-trim(flank(downstream,width=flankSize,start=F))
#
#
#   sm<-ScoreMatrixBin(target=smcRNAseq,windows=anchors,
#                      bin.num=100,
#                      strand.aware=F,weight.col="score",
#                      type="bigWig")
#   smup<-ScoreMatrix(target=smcRNAseq,windows=upstream,
#                     strand.aware=F,weight.col="score",
#                     type="bigwig")
#
#   smdown<-ScoreMatrix(target=smcRNAseq,windows=downstream,
#                       strand.aware=F,weight.col="score",
#                       type="bigwig")
#
#
#   smup2<-ScoreMatrix(target=smcRNAseq,windows=up2,
#                      strand.aware=F,weight.col="score",
#                      type="bigwig")
#
#   smdown2<-ScoreMatrix(target=smcRNAseq,windows=down2,
#                        strand.aware=F,weight.col="score",
#                        type="bigwig")
#
#
#   sml<-as(list(smup2,smup,sm,smdown,smdown2),'ScoreMatrixList')
#   sml<-intersectScoreMatrixList(sml,reorder=T)
#   roworder=rev(order(rowSums(sml[[3]],na.rm=T)))
#   sml<-orderBy(sml,roworder)
#
#   multiHeatMatrix(sml,common.scale=T,winsorize=c(5,95),
#                   matrix.main=c("-20kb","-10kb","anchor","+10kb",
#                                 "+20kb"),
#                   #col=c("blue","lightgrey","red"),
#                   legend.name=prettyGeneName(grp),grid=F)
#
#
#
#   # X chr anchors
#
#   anchors<-keepAnchors
#   anchors<-anchors[seqnames(anchors)=="chrX"]
#   upstream<-trim(flank(anchors,width=flankSize))
#   downstream<-trim(flank(anchors,width=flankSize,start=F))
#
#   upstream<-upstream[order(upstream$region)]
#   downstream<-downstream[order(downstream$region)]
#
#   up2<-trim(flank(upstream,width=flankSize))
#   down2<-trim(flank(downstream,width=flankSize,start=F))
#
#
#   sm<-ScoreMatrixBin(target=smcRNAseq,windows=anchors,
#                      bin.num=100,
#                      strand.aware=F,weight.col="score",
#                      type="bigWig")
#   smup<-ScoreMatrix(target=smcRNAseq,windows=upstream,
#                     strand.aware=F,weight.col="score",
#                     type="bigwig")
#
#   smdown<-ScoreMatrix(target=smcRNAseq,windows=downstream,
#                       strand.aware=F,weight.col="score",
#                       type="bigwig")
#
#
#   smup2<-ScoreMatrix(target=smcRNAseq,windows=up2,
#                      strand.aware=F,weight.col="score",
#                      type="bigwig")
#
#   smdown2<-ScoreMatrix(target=smcRNAseq,windows=down2,
#                        strand.aware=F,weight.col="score",
#                        type="bigwig")
#
#
#   sml<-as(list(smup2,smup,sm,smdown,smdown2),'ScoreMatrixList')
#   sml<-intersectScoreMatrixList(sml,reorder=T)
#   roworder=rev(order(rowSums(sml[[3]],na.rm=T)))
#   sml<-orderBy(sml,roworder)
#
#   multiHeatMatrix(sml,common.scale=T,winsorize=c(5,95),
#                   matrix.main=c("-20kb","-10kb","anchor","+10kb",
#                                 "+20kb"),
#                   #col=c("blue","lightgrey","red"),
#                   legend.name=prettyGeneName(grp),grid=F)
#
#
#
#   # autosomal anchors
#   anchors<-keepAnchors
#   anchors<-anchors[seqnames(anchors)!="chrX"]
#   upstream<-trim(flank(anchors,width=flankSize))
#   downstream<-trim(flank(anchors,width=flankSize,start=F))
#
#   upstream<-upstream[order(upstream$region)]
#   downstream<-downstream[order(downstream$region)]
#
#   up2<-trim(flank(upstream,width=flankSize))
#   down2<-trim(flank(downstream,width=flankSize,start=F))
#
#
#   sm<-ScoreMatrixBin(target=smcRNAseq,windows=anchors,
#                      bin.num=100,
#                      strand.aware=F,weight.col="score",
#                      type="bigWig")
#   smup<-ScoreMatrix(target=smcRNAseq,windows=upstream,
#                     strand.aware=F,weight.col="score",
#                     type="bigwig")
#
#   smdown<-ScoreMatrix(target=smcRNAseq,windows=downstream,
#                       strand.aware=F,weight.col="score",
#                       type="bigwig")
#
#
#   smup2<-ScoreMatrix(target=smcRNAseq,windows=up2,
#                      strand.aware=F,weight.col="score",
#                      type="bigwig")
#
#   smdown2<-ScoreMatrix(target=smcRNAseq,windows=down2,
#                        strand.aware=F,weight.col="score",
#                        type="bigwig")
#
#
#   sml<-as(list(smup2,smup,sm,smdown,smdown2),'ScoreMatrixList')
#   sml<-intersectScoreMatrixList(sml,reorder=T)
#   roworder=rev(order(rowSums(sml[[3]],na.rm=T)))
#   sml<-orderBy(sml,roworder)
#   orderBy(sml,roworder)
#   multiHeatMatrix(sml,common.scale=T,winsorize=c(5,95),
#                   matrix.main=c("-20kb","-10kb","anchor","+10kb",
#                                 "+20kb",
#                                 ylim=c(-2,2)),
#                   #col=c("blue","lightgrey","red"),
#                   legend.name=prettyGeneName(grp),grid=F)
#   dev.off()
# }
#







# anchors - seqplots heatmaps ---------------------------------------------
####
## anchors
####

loops<-import(paste0(outPath,"/otherData/N2.allValidPairs.hic.5-10kbLoops.bedpe"),format="bedpe")

anchors<-zipup(loops)
anchors<-unlist(anchors)
strand(anchors)[seq(1,length(anchors),by=2)]<-"+"
strand(anchors)[seq(2,length(anchors),by=2)]<-"-"

anchors<-reduce(sort(anchors))
#anchors$region<-1:length(anchors)
anchors$chr<-seqnames(anchors)

loopsAll<-paste0(outPath,"/tracks/loops.bed")
export(anchors,con=loopsAll,format="bed")

loopsX<-paste0(outPath,"/tracks/loopsX.bed")
export(anchors[seqnames(anchors)=="chrX"], con=loopsX,format="bed")

loopsA<-paste0(outPath,"/tracks/loopsA.bed")
export(anchors[seqnames(anchors)!="chrX"], con=loopsA,format="bed")

flankSize<-60000

smcRNAseq<-paste0(outPath,"/tracks/salmon_",groupsOI,
                           "_wt_lfc.bw")
if(plotPDFs==T){
  pdf(filename=paste0(outPath,"/plots/anchors-all_flank",flankSize/10000,"kb.pdf"),width=19,
      height=16,units="cm", paper="a4")
}

if(plotPDFs==F){
  png(filename=paste0(outPath,"/plots/anchors-all_flank",flankSize/10000,"kb.png"),width=19,
    height=16,units="cm", res=150)
}

p<-getPlotSetArray(tracks=c(smcRNAseq),
                   features=c(loopsAll),
                refgenome="ce11", bin=100L, xmin=flankSize,
                xmax=flankSize, type="af",
                xanchored=10000)


dd<-plotHeatmap(p,plotz=F)
heatmapQuantiles<-sapply(dd$HLST,quantile,c(0.05,0.95),na.rm=T)
roworder<-rev(order(lapply(dd$HLST,rowSums,na.rm=T)$X1))
minVal<-min(heatmapQuantiles[1,])
maxVal<-max(heatmapQuantiles[2,])
plotHeatmap(p,main="All loop anchors", plotScale="no", sortrows=T,
            clusters=1L,autoscale=F,zmin=minVal, zmax=maxVal,
            indi=F, sort_mids=T,sort_by=c(T,F,F))
if(plotPDFs==F){
  dev.off()
}

if(plotPDFs==F){
  png(filename=paste0(outPath,"/plots/anchors-chrX_flank",flankSize/10000,"kb.png"),width=19,
    height=16,units="cm", res=150)
}
p<-getPlotSetArray(tracks=c(smcRNAseq),
                   features=c(loopsX),
                   refgenome="ce11", bin=100L, xmin=flankSize,
                   xmax=flankSize, type="af",
                   xanchored=10000)
plotHeatmap(p,main="chrX loop anchors", plotScale="no", sortrows=T,
            clusters=1L,autoscale=F,zmin=minVal, zmax=maxVal,
            indi=F, sort_mids=T,sort_by=c(T,F,F))
if(plotPDFs==F){
  dev.off()
}

if(plotPDFs==F){
  png(filename=paste0(outPath,"/plots/anchors-autosomal_flank",flankSize/10000,"kb.png"),width=19,
    height=16,units="cm", res=150)
}

p<-getPlotSetArray(tracks=c(smcRNAseq),
                   features=c(loopsA),
                   refgenome="ce11", bin=100L, xmin=flankSize,
                   xmax=flankSize, type="af",
                   xanchored=10000)
plotHeatmap(p,main="Autosomal loop anchors", plotScale="no", sortrows=T,
            clusters=1L,autoscale=F,zmin=minVal, zmax=maxVal,
            indi=F, sort_mids=T,sort_by=c(T,F,F))
if(plotPDFs==F){
  dev.off()
}

if(plotPDFs==T){
  dev.off()
}

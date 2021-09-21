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
                                    paste0("p",padjVal,"_lfc",lfcVal,"/",
                                           scriptName)))

# AB compartments - N2 ----------------------------------------------------

####
## AB compartments
####

####
## N2 compartments
####
pca2<-import.bw(paste0(outPath,"/otherData/N2_5000b_laminDamID_pca2.bw"))

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
for (grp in useContrasts){
  #grp=useContrasts[1]
  salmon<-readRDS(file=paste0(paste0(outPath,"/rds/",fileNamePrefix,
                  contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))

  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

  salmongr<-sort(salmongr)

  salmongr<-assignGRtoAB(salmongr,pca2,grName=grp,pcaName="N2")
  listgr[[grp]]<-salmongr
}


# plot of counts of significantly changing genes in each compartment
pdf(file=paste0(paste0(outPath,"/plots/",outputNamePrefix,"ABcomp_N2_geneCount_padj",
                       padjVal,"_lfc", lfcVal,".pdf")),
    width=12, height=9, paper="a4r")
par(mfrow=c(2,2))
# genes that change significantly
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="both")

compartmentTable<-do.call(rbind,lapply(lapply(sigList, "[", ,"compartment"),table))
bgCount<-do.call(rbind,lapply(lapply(lapply(listgr,as.data.frame), "[", ,"compartment"),table))

yminmax=c(0,max(compartmentTable))
xx<-barplot(t(compartmentTable),beside=T,col=c("grey80","grey20"),
        main="Significantly changed genes by N2 compartment",cex.axis=1.2,
        cex.names=1.5, ylim=c(yminmax)*1.1)
legend("top",legend = colnames(compartmentTable),fill=c("grey80","grey20"))
text(x=xx, y=t(compartmentTable), label=t(compartmentTable), pos=3,cex=1.3,col="black")


xx<-barplot(t(compartmentTable/bgCount),beside=T,col=c("grey80","grey20"),
            main="Fraction changed genes by N2 compartment",cex.axis=1.2,
            cex.names=1.5,ylim=c(0,0.3))
legend("top",legend = colnames(compartmentTable),fill=c("grey80","grey20"))
text(x=xx, y=t(compartmentTable), label=t(compartmentTable), pos=3,cex=1.3,col="black")

# upregulated genes
sigListUp<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal, lfc=lfcVal, direction="gt")

compartmentTableUp<-do.call(rbind,lapply(lapply(sigListUp, "[", ,"compartment"),table))
colnames(compartmentTableUp)<-paste0(colnames(compartmentTableUp),"_up")


# downregulated genes
sigListDown<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal, lfc= -lfcVal, direction="lt")

compartmentTableDown<-do.call(rbind,lapply(lapply(sigListDown, "[", ,"compartment"),table))
colnames(compartmentTableDown)<-paste0(colnames(compartmentTableDown),"_down")

compartmentTable<-cbind(compartmentTableUp,compartmentTableDown)


Acomp<-compartmentTable[,grep("A",colnames(compartmentTable))]
xx<-barplot(t(Acomp),beside=T,col=c("grey80","grey20"),
        main="Number of up/down regulated in A compartment",cex.axis=1.2,
        cex.names=1.5, ylim=c(yminmax)*1.1)
legend("top", legend=gsub("A_","",colnames(Acomp)), fill=c("grey80","grey20"))
text(x=xx, y=t(Acomp), label=t(Acomp), pos=3,cex=1.3,col="black")


Bcomp<-compartmentTable[,grep("B",colnames(compartmentTable))]
barplot(t(Bcomp),beside=T,col=c("grey80","grey20"),
        main="Number of up/down regulated in B compartment",cex.axis=1.2,
        cex.names=1.5, ylim=c(yminmax)*1.1)
legend("top", legend=gsub("B_","",colnames(Bcomp)), fill=c("grey80","grey20"))
text(x=xx, y=t(Bcomp), label=t(Bcomp), pos=3,cex=1.3,col="black")

dev.off()


####-
## AB compartment by chromosome-----
####-
# genes that change significantly
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="both")

# count genes by category (chr & A/B)
dfl<-lapply(sigList, function(x){x%>% group_by(seqnames,compartment) %>% tally()})
bgCount<-lapply(lapply(listgr,as.data.frame), function(x){x%>% group_by(seqnames,compartment) %>% tally()})
# add name of SMC protein
dfl<-do.call(rbind, mapply(cbind,dfl,"SMC"=names(dfl),SIMPLIFY=F))
bgCount<-do.call(rbind, mapply(cbind,bgCount,"SMC"=names(bgCount),SIMPLIFY=F))

dfl$seqnames<-gsub("chr","",dfl$seqnames)
bgCount$seqnames<-gsub("chr","",bgCount$seqnames)

# do left join to make sure dfl has all the categories required
dfl<-left_join(bgCount,dfl,by=c("seqnames","compartment","SMC"),suffix=c("_total",""))
dfl$n[is.na(dfl$n)]<-0

ymax=max(dfl$n)
p1<-ggplot(dfl,aes(x=seqnames,y=n,group=compartment)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=compartment)) +
  facet_grid(cols=vars(SMC)) +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Significantly changed genes per chromosome by N2 compartment")


dfl$Frac<-dfl$n/dfl$n_total
p1a<-ggplot(dfl,aes(x=seqnames,y=Frac,group=compartment)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=compartment)) +
  facet_grid(cols=vars(SMC)) +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Fraction changed genes per chromosome by N2 compartment")


# upregulated genes
sigListUp<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="gt")
# count genes by category (chr & A/B)
dflUp<-lapply(sigListUp, function(x){x%>% group_by(seqnames,compartment) %>% tally()})
# add name of SMC protein
dflUp<-do.call(rbind, mapply(cbind,dflUp,"SMC"=names(dflUp),SIMPLIFY=F))
dflUp$seqnames<-gsub("chr","",dflUp$seqnames)

# do left join to make sure dfl has all the categories required
dflUp<-left_join(bgCount,dflUp,by=c("seqnames","compartment","SMC"),suffix=c("_total",""))
dflUp$n[is.na(dflUp$n)]<-0
dflUp$expression<-"up"



# downregulated genes
sigListDown<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc= -lfcVal, direction="lt")
dflDown<-lapply(sigListDown, function(x){x%>% group_by(seqnames,compartment) %>% tally()})
# add name of SMC protein
dflDown<-do.call(rbind, mapply(cbind,dflDown,"SMC"=names(dflDown),SIMPLIFY=F))
dflDown$seqnames<-gsub("chr","",dflDown$seqnames)

# do left join to make sure dfl has all the categories required
dflDown<-left_join(bgCount,dflDown,by=c("seqnames","compartment","SMC"),suffix=c("_total",""))
dflDown$n[is.na(dflDown$n)]<-0
dflDown$expression<-"down"

dfl<-rbind(dflUp,dflDown)
dfl$expression<-factor(dfl$expression,levels=c("up","down"))

yminmax=c(0,max(dfl$n))
p2<-ggplot(dfl[dfl$compartment=="A",],aes(x=seqnames,y=n,group=expression)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=expression)) +
  facet_grid(cols=vars(SMC)) +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Up/down regulated in A compartment (N2 pca)")


yminmax=c(0,max(dfl$n))
p3<-ggplot(dfl[dfl$compartment=="B",],aes(x=seqnames,y=n,group=expression)) +
  geom_bar(stat="identity", position=position_dodge(),aes(fill=expression)) +
  facet_grid(cols=vars(SMC),switch="x") +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
  xlab("chr")+ylab("Number of genes") +
  ggtitle("Up/down regulated in B compartment (N2 pca)")


p<-ggpubr::ggarrange(p1,p1a,p2,p3,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                "ABcomp_N2_countsPerChr_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=29,height=19, units="cm")




####-
## AB comp LFC-----
####-

# genes that change significantly
sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="both")
#sigList<-lapply(listgr, as.data.frame)

sigList<-lapply(listgr,as.data.frame)

sigList<-lapply(sigList, "[", ,c("compartment","log2FoldChange"))
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
  ggtitle("Significantly changed genes by N2 compartment") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide="none")

yminmax=c(0,median(abs(sigTbl$log2FoldChange))+quantile(abs(sigTbl$log2FoldChange))[4]*2)
p3<-ggplot(sigTbl,aes(x=compartment,y=abs(log2FoldChange),col=updown,fill=updown)) +
  geom_boxplot(notch=T, varwidth=T, position=position_dodge2(padding=0.2), outlier.shape=NA,
               outlier.color="grey50") +
  facet_grid(cols=vars(SMC)) +
  ggtitle("Significantly changed genes by N2 compartment") +
  theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
  scale_y_continuous(limits = yminmax) +
  scale_color_grey(start=0.2,end=0.2,guide="none")

# test significance of LFC
summary(aov(abs(log2FoldChange)~updown+compartment,data=sigTbl))

p<-ggpubr::ggarrange(p2,p3,ncol=2,nrow=1)
ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                "ABcomp_N2_LFC_padj",
                          padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=29,height=16,units="cm")






# AB compartments - sample specific ---------------------------------------

####
## sample specific compartments -----
####

if(all(c("wt","dpy26cs","kle2cs","scc1cs") %in% varOIlevels)){
  pcas<-data.frame(SMC=c("wt","dpy26cs","kle2cs","scc1cs"),
                   file=list.files(paste0(outPath,"/otherData"),
                                   pattern="_5000_laminDamID_pca2\\.bw"))
  listgr<-NULL
  for (grp in useContrasts){
    #grp=useContrasts[1]
    salmon<-readRDS(file=paste0(outPath,"/rds/",fileNamePrefix,
                contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
    pca2<-import.bw(paste0(outPath,"/otherData/",pcas$file[pcas$SMC==grp]))

    salmon<-salmon[!is.na(salmon$chr),]
    salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

    salmongr<-sort(salmongr)

    salmongr<-assignGRtoAB(salmongr,pca2,grName=grp,pcaName=grp)
    listgr[[grp]]<-salmongr
  }



  pdf(file=paste0(paste0(outPath,"/plots/",outputNamePrefix,
                         "ABcomp_geneCount_padj",
                         padjVal,"_lfc", lfcVal,".pdf")),
      width=19, height=4, paper="a4r")
  par(mfrow=c(1,3))
  # genes that change significantly
  sigList<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                  padj=padjVal,lfc=lfcVal,direction="both")

  compartmentTable<-do.call(rbind,lapply(lapply(sigList, "[", ,"compartment"),table))

  yminmax=c(0,max(compartmentTable))
  xx<-barplot(t(compartmentTable),beside=T,col=c("grey80","grey20"),
              main="Significantly changed genes by compartment",cex.axis=1.2,
              cex.names=1.5, ylim=yminmax*1.1)
  legend("top",legend = colnames(compartmentTable),fill=c("grey80","grey20"))
  text(x=xx, y=t(compartmentTable), label=t(compartmentTable), pos=3,cex=1.3,col="black")


  # upregulated genes
  sigListUp<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                    padj=padjVal,lfc=lfcVal,direction="gt")

  compartmentTableUp<-do.call(rbind,lapply(lapply(sigListUp, "[", ,"compartment"),table))
  colnames(compartmentTableUp)<-paste0(colnames(compartmentTableUp),"_up")


  # downregulated genes
  sigListDown<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                      padj=padjVal, lfc= -lfcVal, direction="lt")

  compartmentTableDown<-do.call(rbind,lapply(lapply(sigListDown, "[", ,"compartment"),table))
  colnames(compartmentTableDown)<-paste0(colnames(compartmentTableDown),"_down")

  compartmentTable<-cbind(compartmentTableUp,compartmentTableDown)



  Acomp<-compartmentTable[,grep("A",colnames(compartmentTable))]
  xx<-barplot(t(Acomp),beside=T,col=c("grey80","grey20"),
              main="Number of up/down regulated in A compartment",cex.axis=1.2,
              cex.names=1.5, ylim=c(yminmax)*1.1)
  legend("top", legend=gsub("A_","",colnames(Acomp)), fill=c("grey80","grey20"))
  text(x=xx, y=t(Acomp), label=t(Acomp), pos=3,cex=1.3,col="black")


  Bcomp<-compartmentTable[,grep("B",colnames(compartmentTable))]
  barplot(t(Bcomp),beside=T,col=c("grey80","grey20"),
          main="Number of up/down regulated in B compartment",cex.axis=1.2,
          cex.names=1.5, ylim=c(yminmax)*1.1)
  legend("top", legend=gsub("B_","",colnames(Bcomp)), fill=c("grey80","grey20"))
  text(x=xx, y=t(Bcomp), label=t(Bcomp), pos=3,cex=1.3,col="black")

  dev.off()


  ####
  ## AB compartment by chromosome -----
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
  sigListUp<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                    padj=padjVal,lfc=lfcVal,direction="gt")
  # count genes by category (chr & A/B)
  dflUp<-lapply(sigListUp, function(x){x%>% group_by(seqnames,compartment) %>% tally()})
  # add name of SMC protein
  dflUp<-do.call(rbind, mapply(cbind,dflUp,"SMC"=names(dflUp),SIMPLIFY=F))
  dflUp$seqnames<-gsub("chr","",dflUp$seqnames)
  dflUp$expression<-"up"

  # downregulated genes
  sigListDown<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                      padj=padjVal, lfc= -lfcVal, direction="lt")
  dflDown<-lapply(sigListDown, function(x){x%>% group_by(seqnames,compartment) %>% tally()})
  # add name of SMC protein
  dflDown<-do.call(rbind, mapply(cbind,dflDown,"SMC"=names(dflDown),SIMPLIFY=F))
  dflDown$seqnames<-gsub("chr","",dflDown$seqnames)
  dflDown$expression<-"down"

  dfl<-rbind(dflUp,dflDown)
  dfl$expression<-factor(dfl$expression,levels=c("up","down"))

  yminmax=c(0,max(dfl$n))
  p2<-ggplot(dfl[dfl$compartment=="A",],aes(x=seqnames,y=n,group=expression)) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=expression)) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
    xlab("chr")+ylab("Number of genes") +
    ggtitle("Up/down regulated in A compartment")


  yminmax=c(0,max(dfl$n))
  p3<-ggplot(dfl[dfl$compartment=="B",],aes(x=seqnames,y=n,group=expression)) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=expression)) +
    facet_grid(cols=vars(SMC),switch="x") +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
    xlab("chr")+ylab("Number of genes") +
    ggtitle("Up/down regulated in B compartment")


  p<-ggpubr::ggarrange(p1,p2,p3,ncol=1,nrow=3)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                  "ABcomp_countsPerChr_padj",
                                  padjVal,"_lfc", lfcVal,".pdf"),
                  plot=p, device="pdf",width=19,height=29,units="cm")



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
    scale_y_continuous(limits = yminmax) +
    scale_color_grey(start=0.2,end=0.2,guide="none")

  yminmax=c(0,median(abs(sigTbl$log2FoldChange))+quantile(abs(sigTbl$log2FoldChange))[4]*2)
  p3<-ggplot(sigTbl,aes(x=compartment,y=abs(log2FoldChange),col=updown,fill=updown)) +
    geom_boxplot(notch=T, varwidth=T, position=position_dodge2(padding=0.2),
                 outlier.shape=NA) +
    facet_grid(cols=vars(SMC)) +
    ggtitle("Significantly changed genes by compartment") +
    theme_minimal() + scale_fill_grey(start=0.8,end=0.3) +
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

  p<-ggpubr::ggarrange(p1,p1a,ncol=1,nrow=3)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                  "ABcompSwitch_countsPerChr_ABBA_padj",
                                  padjVal,"_lfc", lfcVal,".pdf"),
                  plot=p, device="pdf",width=19,height=29,units="cm")


  # upregulated genes
  sigListUp<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                    padj=padjVal,lfc=lfcVal,direction="gt")
  # count genes by category (chr & A/B)
  dflUp<-lapply(sigListUp, function(x){x%>% group_by(seqnames,switch,.drop=F) %>% tally()})
  # add name of SMC protein
  dflUp<-do.call(rbind, mapply(cbind,dflUp,"SMC"=names(dflUp),SIMPLIFY=F))
  dflUp$seqnames<-gsub("chr","",dflUp$seqnames)
  dflUp$expression<-"up"

  # downregulated genes
  sigListDown<-lapply(lapply(listgr,as.data.frame), getSignificantGenes,
                      padj=padjVal, lfc= -lfcVal, direction="lt")
  dflDown<-lapply(sigListDown, function(x){x%>% group_by(seqnames,switch,.drop=F) %>% tally()})
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
    xlab("chr")+ylab("Number of genes") + ylim(yminmax) +
    ggtitle("Up/down regulated genes in AA compartment")

  p3<-ggplot(dfl[dfl$switch=="BB",],aes(x=seqnames,y=n,group=expression)) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=expression)) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
    xlab("chr")+ylab("Number of genes") + ylim(yminmax) +
    ggtitle("Up/down regulated genes in BB compartment")

  yminmax=c(0,max(dfl$n[! (dfl$switch %in% c("AA","BB"))]))
  p4<-ggplot(dfl[dfl$switch=="AB",],aes(x=seqnames,y=n,group=expression)) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=expression)) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
    xlab("chr")+ylab("Number of genes") + ylim(yminmax) +
    ggtitle("Up/down regulated genes in AB compartment")

  p5<-ggplot(dfl[dfl$switch=="BA",],aes(x=seqnames,y=n,group=expression)) +
    geom_bar(stat="identity", position=position_dodge(),aes(fill=expression)) +
    facet_grid(cols=vars(SMC)) +
    theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
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
    scale_y_continuous(limits = yminmax) +
    scale_color_grey(start=0.7,end=0.3,guide="none")



  p<-ggpubr::ggarrange(p1,p2,ncol=2,nrow=1)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
                                  "ABcompSwitch_LFC_padj",
                                  padjVal,"_lfc", lfcVal,".pdf"),
                  plot=p, device="pdf",width=29,height=16,units="cm")


}


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
# for (grp in useContrasts){
#   anchors<-keepAnchors
#   smcRNAseq<-import(paste0(outPath,"/tracks/",outputNamePrefix,grp,
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

loops<-rtracklayer::import(paste0(outPath,"/otherData/N2.allValidPairs.hic.5-10kbLoops.bedpe"),
                           format="bedpe")

anchors<-zipup(loops)
anchors<-unlist(anchors)
strand(anchors)[seq(1,length(anchors),by=2)]<-"+"
strand(anchors)[seq(2,length(anchors),by=2)]<-"-"

anchors<-reduce(sort(anchors))
#anchors$region<-1:length(anchors)
anchors$chr<-seqnames(anchors)

loopsAll<-paste0(outPath,"/tracks/",outputNamePrefix,"loops.bed")
export(anchors,con=loopsAll,format="bed")

loopsX<-paste0(outPath,"/tracks/",outputNamePrefix,"loopsX.bed")
export(anchors[seqnames(anchors)=="chrX"], con=loopsX,format="bed")

loopsA<-paste0(outPath,"/tracks/",outputNamePrefix,"loopsA.bed")
export(anchors[seqnames(anchors)!="chrX"], con=loopsA,format="bed")

flankSize<-60000

smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
                  useContrasts,"_lfc.bw")
if(plotPDFs==T){
  pdf(file=paste0(outPath,"/plots/",outputNamePrefix,"anchors-all_flank",
                      flankSize/1000,"kb.pdf"), width=19,
      height=16, paper="a4")
}

if(plotPDFs==F){
  png(filename=paste0(outPath,"/plots/",outputNamePrefix,"anchors-all_flank",
                      flankSize/1000,"kb.png"),width=19,
    height=16, units="cm", res=150)
}

#############################
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
  png(filename=paste0(outPath,"/plots/",outputNamePrefix,"anchors-chrX_flank",
                      flankSize/1000,"kb.png"), width=19,
    height=16, units="cm", res=150)
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
  png(filename=paste0(outPath,"/plots/",outputNamePrefix,"anchors-autosomal_flank",
                      flankSize/1000,"kb.png"), width=19,
    height=16, units="cm", res=150)
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


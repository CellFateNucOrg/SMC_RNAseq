library(readxl)
library(ggVennDiagram)
library(ggplot2)
library(EnhancedVolcano)

source("functions.R")

outPath="."
padjVal=0.05
lfcVal=0.5
plotPDFs=T
fileNamePrefix="noOsc_"
filterPrefix="noOsc_"
filterData=T

fileList<-read.table(paste0(outPath,"/fastqList.txt"),stringsAsFactors=F,header=T)

# extract the strain variable
strain<-factor(as.character(unique(fileList$sampleName),levels=c("366","382","775","784")))
SMC<-strain
levels(SMC)<-c("wt","dpy26cs","kle2cs","scc1cs")

controlGrp<-levels(SMC)[1] # control group
groupsOI<-levels(SMC)[-1]


###########################
## compare samples
##########################


#######
## venn diagrams
#######

sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                                                      namePadjCol="padj",
                                                                      nameLfcCol="log2FoldChange",
                                                                      direction="both",
                                                                      chr="all", nameChrCol="chr"))
}

sigGenes<-lapply(sigTables, "[", ,"wormbaseID")




for (grp in groupsOI){
  #grp="dpy26cs"
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

  if(filterData){
    oscillating<-read.delim(paste0(outPath,"/oscillatingGenes.tsv"),header=T,
                            stringsAsFactors=F) #3739
    #latorre<-read.delim(paste0(outPath,"/oscillatingGenes_latorre.tsv")) #3235
    #hsUP<-readRDS(file="hsUp_garrigues2019.rds") #1680
    #hsDOWN<-readRDS(file="hsDown_garrigues2019.rds") #455

    toFilter<-unique(c(oscillating$WB_ID))
    #toFilter<-unique(c(oscillating$WB_ID,latorre$wormbaseID,hsUP$WormBase.ID,
    #                   hsDOWN$WormBase.ID))
    #4522 genes osc+latorre
    #6101 genes osc+latorre+hs

    # remove filtered genes
    idx<-salmon$wormbaseID %in% toFilter
    salmon<-salmon[!idx,]

    fileNamePrefix=filterPrefix
  }

  #####################################################
  ## compare to public data
  #####################################################

  kramerURL<-"https://doi.org/10.1371/journal.pgen.1005698.s011"
  kramerFileName="Kramer_2015_PlotGen_S3file.xlsx"
  if(!file.exists(kramerFileName)){
    download.file(url=kramerURL,destfile=kramerFileName)
  }
  kramer<-read_excel(kramerFileName,col_types=c(rep("text",3),rep("numeric",30)))
  names(kramer)
  kramer<-kramer[,c(1:3,grep("_L3_",names(kramer)))]


  kramerDpy27<-getSignificantGenes(kramer, padj=padjVal, lfc=lfcVal,
                                   namePadjCol="dpy27_RNAi_L3_padj",
                                   nameLfcCol="dpy27_RNAi_L3_log2_fold_change",
                                   direction="both",
                                   chr="all", nameChrCol="chr", outPath=".")

  kramerDpy21<-getSignificantGenes(kramer, padj=padjVal, lfc=lfcVal,
                                   namePadjCol="dpy21_mutant_L3_padj",
                                   nameLfcCol="dpy21_mutant_L3_log2_fold_change",
                                   direction="both",
                                   chr="all", nameChrCol="chr")

  if(filterData){
    # remove filtered genes
    idx<-kramerDpy27$Gene_WB_ID %in% toFilter
    kramerDpy27<-kramerDpy27[!idx,]

    idx<-kramerDpy21$Gene_WB_ID %in% toFilter
    kramerDpy21<-kramerDpy21[!idx,]
  }

  salmonSig<-getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                 namePadjCol="padj",
                                 nameLfcCol="log2FoldChange",
                                 direction="both",
                                 chr="all", nameChrCol="chr")


  x<-list(salmon=salmonSig$wormbaseID, dpy27=kramerDpy27$Gene_WB_ID,
          dpy21=kramerDpy21$Gene_WB_ID)
  names(x)<-c(prettyGeneName(grp), "dpy-27", "dpy-21")
  txtLabels<-list()
  txtLabels[paste0("% ",names(x)[2]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[2]]))/length(x[[2]]),1)
  txtLabels[paste0("% ",names(x)[3]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[3]]),1)
  p1<-ggVennDiagram(x) + ggtitle(label=paste0(grp," vs Kramer(2015): |lfc|>", lfcVal, ", padj<",padjVal),subtitle=paste0(txtLabels[[1]],names(txtLabels)[1], " & ", txtLabels[[2]], names(txtLabels)[2]))

  #p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_",grp,"VsKarmer_padj",
                                  padjVal,"_lfc", lfcVal,".pdf"),
                  plot=p1, device="pdf",width=15,height=11,units="cm")



  ###############################
  ## chrX upregulated
  ###############################

  dim(kramer)
  dim(salmon)
  idx<-match(kramer$Gene_WB_ID, salmon$wormbaseID)
  kramer$chr<-salmon$chr[idx]
  kramer<-kramer[!is.na(kramer$chr),]


  salmondc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")

  idx<-!is.na(kramer$dpy27_RNAi_L3_padj) &
    kramer$dpy27_RNAi_L3_padj < padjVal &
    kramer$dpy27_RNAi_L3_log2_fold_change > lfcVal &
    kramer$chr=="chrX"
  kramerdpy27dc<-kramer[idx,]


  idx<-!is.na(kramer$dpy21_mutant_L3_padj) &
    kramer$dpy21_mutant_L3_padj < padjVal &
    kramer$dpy21_mutant_L3_log2_fold_change > lfcVal &
    kramer$chr=="chrX"
  kramerdpy21dc<-kramer[idx,]

  x<-list(salmon=salmondc$wormbaseID, dpy27=kramerdpy27dc$Gene_WB_ID,
          dpy21=kramerdpy21dc$Gene_WB_ID)
  names(x)<-c(prettyGeneName(grp), "dpy-27", "dpy-21")
  txtLabels<-list()
  txtLabels[paste0("% ",names(x)[2]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[2]]))/length(x[[2]]),1)
  txtLabels[paste0("% ",names(x)[3]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[3]]),1)

  p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DC genes: lfc>", lfcVal, ", padj<",padjVal), subtitle=paste0(txtLabels[[1]],names(txtLabels)[1], " & ", txtLabels[[2]], names(txtLabels)[2]))


  #p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,
                                  "venn_chrXup_",grp, "VsKramer_padj",
                                  padjVal,"_lfc", lfcVal,".pdf"),
                  plot=p1, device="pdf",width=15,height=11,units="cm")


  ###############################
  ## autosome changed
  ###############################

  dim(kramer)
  dim(salmon)
  idx<-match(kramer$Gene_WB_ID, salmon$wormbaseID)
  kramer$chr<-salmon$chr[idx]
  kramer<-kramer[!is.na(kramer$chr),]


  salmondc<-filterResults(salmon,padjVal,lfcVal,direction="both",chr="autosomes")

  idx<-!is.na(kramer$dpy27_RNAi_L3_padj) &
    kramer$dpy27_RNAi_L3_padj < padjVal &
    abs(kramer$dpy27_RNAi_L3_log2_fold_change) > lfcVal &
    kramer$chr!="chrX"
  kramerdpy27dc<-kramer[idx,]


  idx<-!is.na(kramer$dpy21_mutant_L3_padj) &
    kramer$dpy21_mutant_L3_padj < padjVal &
    abs(kramer$dpy21_mutant_L3_log2_fold_change) > lfcVal &
    kramer$chr!="chrX"
  kramerdpy21dc<-kramer[idx,]

  x<-list(salmon=salmondc$wormbaseID, dpy27=kramerdpy27dc$Gene_WB_ID,
          dpy21=kramerdpy21dc$Gene_WB_ID)
  names(x)<-c(prettyGeneName(grp), "dpy-27", "dpy-21")
  txtLabels<-list()
  txtLabels[paste0("% ",names(x)[2]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[2]]))/length(x[[2]]),1)
  txtLabels[paste0("% ",names(x)[3]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[3]]),1)
  p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("Autosomal genes: |lfc|>", lfcVal, ", padj<",padjVal),subtitle=paste0(txtLabels[[1]],names(txtLabels)[1], " & ", txtLabels[[2]], names(txtLabels)[2]))


  #p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,
                                  "venn_chrAall_",grp, "VsKramer_padj",
                                  padjVal,"_lfc", lfcVal,".pdf"),
                  plot=p1, device="pdf",width=15,height=11,units="cm")

}


grp="dpy26cs"

###############################
## classical dc genes
###############################

pubDC<-readRDS("/Users/semple/Documents/MeisterLab/dSMF/DCgenes/published_DCgr.rds")
pubNDC<-readRDS("/Users/semple/Documents/MeisterLab/dSMF/DCgenes/published_NDCgr.rds")

if(filterData){
  # remove filtered genes
  idx<-pubDC$wormbase %in% toFilter
  pubDC<-pubDC[!idx,]

  idx<-pubNDC$wormbase %in% toFilter
  pubNDC<-pubNDC[!idx,]
}


salmondc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondc$wormbaseID, DC=pubDC$wormbase,
        nonDC=pubNDC$wormbase)
names(x)<-c(prettyGeneName(grp), "DC", "nonDC")
txtLabels<-list()
txtLabels[paste0("% ",names(x)[2]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[2]]))/length(x[[2]]),1)
txtLabels[paste0("% ",names(x)[3]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[3]]),1)
p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX dc genes: lfc>", lfcVal, ", padj<",padjVal), subtitle=paste0(txtLabels[[1]],names(txtLabels)[1], " & ", txtLabels[[2]], names(txtLabels)[2]))


#p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_",
                                grp, "VsPapers_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p1, device="pdf",width=15,height=11,units="cm")






###############################
## Jans 2009
###############################

JansDC<-readRDS("/Users/semple/Documents/MeisterLab/dSMF/DCgenes/Jans2009_DCgr.rds")
JansNDC<-readRDS("/Users/semple/Documents/MeisterLab/dSMF/DCgenes/Jans2009_NDCgr.rds")

if(filterData){
  # remove filtered genes
  idx<-JansDC$wormbase %in% toFilter
  JansDC<-JansDC[!idx,]

  idx<-JansNDC$wormbase %in% toFilter
  JansNDC<-JansNDC[!idx,]
}



salmondc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondc$wormbaseID, JansDC=JansDC$wormbase,
        JansNDC=JansNDC$wormbase)
names(x)<-c(prettyGeneName(grp), "JansDC", "JansNDC")
txtLabels<-list()
txtLabels[paste0("% ",names(x)[2]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[2]]))/length(x[[2]]),1)
txtLabels[paste0("% ",names(x)[3]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[3]]),1)
p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX dc genes: lfc>", lfcVal, ", padj<",padjVal),subtitle=paste0(txtLabels[[1]],names(txtLabels)[1], " & ", txtLabels[[2]], names(txtLabels)[2]))

#p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_",grp,"VsJans2009_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p1, device="pdf",width=15,height=11,units="cm")




###############################
## Jans 2009 vs Kramer
###############################

kramer<-read_excel(kramerFileName,col_types=c(rep("text",3),rep("numeric",30)))
names(kramer)
kramer<-kramer[,c(1:3,grep("_L3_",names(kramer)))]

idx<-!is.na(kramer$dpy27_RNAi_L3_padj) & abs(kramer$dpy27_RNAi_L3_log2_fold_change)>lfcVal & kramer$dpy27_RNAi_L3_padj<padjVal
kramerDpy27<-kramer[idx,]

idx<-!is.na(kramer$dpy21_mutant_L3_padj) & abs(kramer$dpy21_mutant_L3_log2_fold_change)>lfcVal & kramer$dpy21_mutant_L3_padj<padjVal
kramerDpy21<-kramer[idx,]


x<-list(JansDC=JansDC$wormbase, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)
names(x)<-c("JansDC", "dpy-27", "dpy-21")
txtLabels<-list()
txtLabels[paste0("% ",names(x)[2]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[2]]))/length(x[[2]]),1)
txtLabels[paste0("% ",names(x)[3]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[3]]),1)
p1<-ggVennDiagram(x) + ggtitle(label=paste0("Jans DC vs Kramer(2015): lfc>", lfcVal, ", padj<",padjVal), subtitle=paste0(txtLabels[[1]],names(txtLabels)[1], " & ", txtLabels[[2]], names(txtLabels)[2]))


x<-list(JansNDC=JansNDC$wormbase, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)
names(x)<-c("JansNDC", "dpy-27", "dpy-21")
txtLabels<-list()
txtLabels[paste0("% ",names(x)[2]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[2]]))/length(x[[2]]),1)
txtLabels[paste0("% ",names(x)[3]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[3]]),1)
p2<-ggVennDiagram(x) + ggtitle(label=paste0("Jans NDC vs Kramer(2015): lfc>", lfcVal, ", padj<",padjVal), subtitle=paste0(txtLabels[[1]],names(txtLabels)[1], " & ", txtLabels[[2]], names(txtLabels)[2]))



p<-ggpubr::ggarrange(p1,p2,ncol=2,nrow=1)
ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,
                                "venn_Jans2009vKramer2015_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=29,height=11,units="cm")
<<<<<<< HEAD



###############-
## intersection oscillating -----
###############-
if(!filterData){
  oscillating<-read.delim(paste0(outPath,"/oscillatingGenes.tsv"),header=T,
                          stringsAsFactors=F) #3739
  latorre<-read.delim(paste0(outPath,"/oscillatingGenes_latorre.tsv")) #3235
  hsUP<-readRDS(file="hsUp_garrigues2019.rds") #1680
  hsDOWN<-readRDS(file="hsDown_garrigues2019.rds") #455

  for (grp in groupsOI){
    salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))
    salmonSig<-getSignificantGenes(salmon, padj=padjVal,
                                   lfc=lfcVal,
                                   namePadjCol="padj",
                                   nameLfcCol="log2FoldChange",
                                   direction="both",
                                   chr="all", nameChrCol="chr")
=======
>>>>>>> a45e54b (add subtitle with set overlap)

    x<-list(salmon=salmonSig$wormbaseID, Meeuse=oscillating$WB_ID,
            Latorre=latorre$wormbaseID)
    names(x)<-c(prettyGeneName(grp), "Meeuse", "Latorre")
    txtLabels<-list()
    txtLabels[paste0("% ",names(x)[1]," in ",names(x)[2])]<-round(100*length(intersect(x[[1]],x[[2]]))/length(x[[1]]),1)
    txtLabels[paste0("% ",names(x)[1]," in ",names(x)[3])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[1]]),1)
    txtLabels[paste0("% ",names(x)[1]," in union")]<-round(100*length(intersect(x[[1]],union(x[[2]],x[[3]])))/length(x[[1]]),1)
    p1<-ggVennDiagram(x) + ggtitle(label=paste0(grp," vs oscillating: |lfc|>", lfcVal, ", padj<",padjVal),subtitle=paste0(txtLabels[[1]],names(txtLabels)[1], " & ", txtLabels[[2]], names(txtLabels)[2], "\n",txtLabels[[3]],names(txtLabels)[3]))

    #p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
    ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_",grp,"VsOscillating_padj",
                                    padjVal,"_lfc", lfcVal,".pdf"),
                    plot=p1, device="pdf",width=15,height=11,units="cm")

    # hsUP
    salmonUp<-getSignificantGenes(salmon, padj=padjVal,
                                   lfc=lfcVal,
                                   namePadjCol="padj",
                                   nameLfcCol="log2FoldChange",
                                   direction="gt",
                                   chr="all", nameChrCol="chr")
    salmonDown<-getSignificantGenes(salmon, padj=padjVal,
                                  lfc=lfcVal,
                                  namePadjCol="padj",
                                  nameLfcCol="log2FoldChange",
                                  direction="lt",
                                  chr="all", nameChrCol="chr")

    x<-list(up=salmonUp$wormbaseID, down=salmonDown$wormbaseID,
            hs_up=hsUP$WormBase.ID, hs_down=hsDOWN$WormBase.ID)
    names(x)<-c(paste0(prettyGeneName(grp),c("_up","_down")),
                       "hs_up", "hs_down")
    txtLabels<-list()
    txtLabels[paste0("% ",names(x)[1]," in ",names(x)[3])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[1]]),1)
    txtLabels[paste0("% ",names(x)[1]," in ",names(x)[4])]<-round(100*length(intersect(x[[1]],x[[4]]))/length(x[[1]]),1)
    txtLabels[paste0("% ",names(x)[2]," in ",names(x)[3])]<-round(100*length(intersect(x[[2]],x[[3]]))/length(x[[2]]),1)
    txtLabels[paste0("% ",names(x)[2]," in ",names(x)[4])]<-round(100*length(intersect(x[[2]],x[[4]]))/length(x[[2]]),1)
    p1<-ggVennDiagram(x) + ggtitle(label=paste0(grp," vs oscillating: |lfc|>", lfcVal, ", padj<",padjVal),subtitle=paste0(txtLabels[[1]], names(txtLabels)[1], " & ", txtLabels[[2]], names(txtLabels)[2], " & \n", txtLabels[[3]], names(txtLabels)[3], " & ", txtLabels[[4]], names(txtLabels)[4]))

    #p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
    ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_",grp,"VsHeatshock_padj",
                                    padjVal,"_lfc", lfcVal,".pdf"),
                    plot=p1, device="pdf",width=15,height=11,units="cm")

  }
}

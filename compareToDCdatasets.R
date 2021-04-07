library(readxl)
library(ggVennDiagram)
library(ggplot2)
library(EnhancedVolcano)

source("functions.R")
source("./variableSettings.R")
if(filterData){
  fileNamePrefix=filterPrefix
}

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
sigGR<-list()
df<-data.frame(matrix(nrow=10,ncol=3))
names(df)<-groupsOI
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))

  n<-prettyGeneName(grp)
  sigTables[[n]]<-as.data.frame(getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                                    namePadjCol="padj",
                                                    nameLfcCol="log2FoldChange",
                                                    direction="both",
                                                    chr="all", nameChrCol="chr"))
  sigTables[[n]]<-sigTables[[n]][!is.na(sigTables[[n]]$chr),] # removes mtDNA genes
  sigGR[[n]]<-GRanges(seqnames=sigTables[[n]]$chr,
                      ranges=IRanges(start=sigTables[[n]]$start,
                                     end=sigTables[[n]]$end),
                      strand=sigTables[[n]]$strand)
  mcols(sigGR[[n]])<-sigTables[[n]]

  # check if datasets have chrX genes included
  includeChrX<-"chrX" %in% unlist(lapply(sigTables,"[","chr"))
  if(includeChrX){
    chrXgr<-sigGR[[n]][sigGR[[n]]$chr=="chrX"]
    print(chrXgr[order(width(chrXgr),decreasing=T)][1:10])
    print(chrXgr[order(width(chrXgr),decreasing=T)]$wormbaseID[1:10])
    df[,grp]<-print(chrXgr[order(width(chrXgr),decreasing=T)]$publicID[1:10])
    print(width(chrXgr[order(width(chrXgr),decreasing=T)][1:10]))
    hist(width(chrXgr),main=n,breaks=100)
  }
}
df


for (grp in groupsOI){
  #grp="dpy26cs"
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))

  if(filterData){
    # remove filtered genes
    idx<-salmon$wormbaseID %in% toFilter
    salmon<-salmon[!idx,]

    fileNamePrefix=filterPrefix
  }

  #####################################################
  ## compare to public data
  #####################################################

  # kramerURL<-"https://doi.org/10.1371/journal.pgen.1005698.s011"
  # kramerFileName="Kramer_2015_PlotGen_S3file.xlsx"
  # if(!file.exists(paste0(outPath,"/publicData/",kramerFileName))){
  #   download.file(url=kramerURL,
  #                 destfile=paste0(outPath,"/publicData/",kramerFileName))
  # }
  # kramer<-read_excel(paste0(outPath,"/publicData/",kramerFileName),
  #                    col_types=c(rep("text",3),rep("numeric",30)))
  # names(kramer)
  # kramer<-kramer[,c(1:3,grep("_L3_",names(kramer)))]
  kramer<-as.data.frame(readRDS(file=paste0(outPath,"/publicData/kramer2015_L3_gr.rds")))

  localPadj=0.05
  localLFC=0.25
  kramerDpy27<-getSignificantGenes(kramer, padj=localPadj, lfc=localLFC,
                                   namePadjCol="dpy27_RNAi_L3_padj",
                                   nameLfcCol="dpy27_RNAi_L3_log2_fold_change",
                                   direction="both",
                                   chr="all", nameChrCol="seqnames", outPath=".")

  kramerDpy21<-getSignificantGenes(kramer, padj=localPadj, lfc=localLFC,
                                   namePadjCol="dpy21_mutant_L3_padj",
                                   nameLfcCol="dpy21_mutant_L3_log2_fold_change",
                                   direction="both",
                                   chr="all", nameChrCol="seqnames")

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


  x<-list(salmon=salmonSig$wormbaseID, dpy27=kramerDpy27$wormbaseID,
          dpy21=kramerDpy21$wormbaseID)
  names(x)<-c(prettyGeneName(grp), "dpy-27", "dpy-21")
  txtLabels<-list()
  txtLabels[paste0("% ",names(x)[2]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[2]]))/length(x[[2]]),1)
  txtLabels[paste0("% ",names(x)[3]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[3]]),1)
  p1<-ggVennDiagram(x) + ggtitle(label=paste0(grp," vs Kramer(2015): |lfc|>", lfcVal, ", padj<",padjVal),subtitle=paste0(txtLabels[[1]],names(txtLabels)[1], " & ", txtLabels[[2]], names(txtLabels)[2]))

  #p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_",
                                  grp,"VsKarmer_padj",
                                  padjVal,"_lfc", lfcVal,".pdf"),
                  plot=p1, device="pdf",width=15,height=11,units="cm")


  if(includeChrX){
    ###############################-
    ## chrX upregulated-----
    ###############################-

    dim(kramer)
    dim(salmon)
    #idx<-match(kramer$wormbaseID, salmon$wormbaseID)
    #kramer$chr<-salmon$chr[idx]
    #kramer<-kramer[!is.na(kramer$chr),]


    salmondc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")

    kramerdpy27dc<-getSignificantGenes(kramer, padj=localPadj, lfc=localLFC,
                                     namePadjCol="dpy27_RNAi_L3_padj",
                                     nameLfcCol="dpy27_RNAi_L3_log2_fold_change",
                                     direction="gt",
                                     chr="chrX", nameChrCol="seqnames",
                                     outPath=outPath)

    kramerdpy21dc<-getSignificantGenes(kramer, padj=localPadj, lfc=localLFC,
                                     namePadjCol="dpy21_mutant_L3_padj",
                                     nameLfcCol="dpy21_mutant_L3_log2_fold_change",
                                     direction="gt",
                                     chr="chrX", nameChrCol="seqnames",
                                     outPath=outPath)

    x<-list(salmon=salmondc$wormbaseID, dpy27=kramerdpy27dc$wormbaseID,
            dpy21=kramerdpy21dc$wormbaseID)
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
  }

  ###############################
  ## autosome changed
  ###############################

  dim(kramer)
  dim(salmon)
  #idx<-match(kramer$wormbaseID, salmon$wormbaseID)
  #kramer$chr<-salmon$chr[idx]
  #kramer<-kramer[!is.na(kramer$chr),]


  salmondc<-filterResults(salmon,padjVal,lfcVal,direction="both",chr="autosomes")

  kramerdpy27dc<-getSignificantGenes(kramer, padj=localPadj, lfc=localLFC,
                                     namePadjCol="dpy27_RNAi_L3_padj",
                                     nameLfcCol="dpy27_RNAi_L3_log2_fold_change",
                                     direction="both",
                                     chr="autosomes", nameChrCol="seqnames",
                                     outPath=outPath)

  kramerdpy21dc<-getSignificantGenes(kramer, padj=localPadj, lfc=localLFC,
                                     namePadjCol="dpy21_mutant_L3_padj",
                                     nameLfcCol="dpy21_mutant_L3_log2_fold_change",
                                     direction="both",
                                     chr="autosomes", nameChrCol="seqnames",
                                     outPath=outPath)

  x<-list(salmon=salmondc$wormbaseID, dpy27=kramerdpy27dc$worbaseID,
          dpy21=kramerdpy21dc$wormbaseID)
  names(x)<-c(prettyGeneName(grp), "dpy-27", "dpy-21")
  txtLabels<-list()
  txtLabels[paste0("% ",names(x)[2]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[2]]))/length(x[[2]]),1)
  txtLabels[paste0("% ",names(x)[3]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[3]]),1)
  p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("Autosomal genes: |lfc|>", lfcVal, ", padj<",padjVal), subtitle=paste0(txtLabels[[1]], names(txtLabels)[1], " & ", txtLabels[[2]], names(txtLabels)[2]))


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

pubDC<-readRDS(paste0(outPath,"/publicData/published_DCgr.rds"))
pubNDC<-readRDS(paste0(outPath,"/publicData/published_NDCgr.rds"))

if(filterData){
  # remove filtered genes
  idx<-pubDC$wormbaseID %in% toFilter
  pubDC<-pubDC[!idx,]

  idx<-pubNDC$wormbaseID %in% toFilter
  pubNDC<-pubNDC[!idx,]
}

if(includeChrX){
  salmondc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
  x<-list(salmon=salmondc$wormbaseID, DC=pubDC$wormbaseID,
          nonDC=pubNDC$wormbaseID)
  names(x)<-c(prettyGeneName(grp), "DC", "nonDC")
  txtLabels<-list()
  txtLabels[paste0("% ",names(x)[2]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[2]]))/length(x[[2]]),1)
  txtLabels[paste0("% ",names(x)[3]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[3]]),1)
  p1<-ggVennDiagram(x) +
    ggplot2::ggtitle(label=paste0("chrX dc genes: lfc>", lfcVal, ", padj<",
                                  padjVal),
                     subtitle=paste0(txtLabels[[1]], names(txtLabels)[1],
                                  " & ", txtLabels[[2]], names(txtLabels)[2]))


  #p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_",
                                  grp, "VsPapers_padj",
                                  padjVal,"_lfc", lfcVal,".pdf"),
                  plot=p1, device="pdf",width=15,height=11,units="cm")
}





###############################
## Jans 2009
###############################

JansDC<-readRDS(paste0(outPath,"/publicData/Jans2009_DCgr.rds"))
JansNDC<-readRDS(paste0(outPath,"/publicData/Jans2009_NDCgr.rds"))

if(filterData){
  # remove filtered genes
  idx<-JansDC$wormbaseID %in% toFilter
  JansDC<-JansDC[!idx,]

  idx<-JansNDC$wormbaseID %in% toFilter
  JansNDC<-JansNDC[!idx,]
}


if(includeChrX){
  salmondc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
  x<-list(salmon=salmondc$wormbaseID, JansDC=JansDC$wormbaseID,
          JansNDC=JansNDC$wormbaseID)
  names(x)<-c(prettyGeneName(grp), "JansDC", "JansNDC")
  txtLabels<-list()
  txtLabels[paste0("% ",names(x)[2]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[2]]))/length(x[[2]]),1)
  txtLabels[paste0("% ",names(x)[3]," in ",names(x)[1])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[3]]),1)
  p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX dc genes: lfc>", lfcVal, ", padj<",padjVal),subtitle=paste0(txtLabels[[1]],names(txtLabels)[1], " & ", txtLabels[[2]], names(txtLabels)[2]))

  #p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_",grp,
                                  "VsJans2009_padj", padjVal,"_lfc", lfcVal,
                                  ".pdf"),
                  plot=p1, device="pdf",width=15,height=11,units="cm")
}





###############-
## intersection oscillating -----
###############-
if(!filterData){
  oscillating<-read.delim(paste0(outPath,"/publicData/oscillatingGenes.tsv"),header=T,
                          stringsAsFactors=F) #3739
  latorre<-read.delim(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv")) #3235
  hsUP<-readRDS(file=paste0(outPath,"/publicData/hsUp_garrigues2019.rds")) #1680
  hsDOWN<-readRDS(file=paste0(outPath,"/publicData/hsDown_garrigues2019.rds")) #455

  for (grp in groupsOI){
    salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))
    salmonSig<-getSignificantGenes(salmon, padj=padjVal,
                                   lfc=lfcVal,
                                   namePadjCol="padj",
                                   nameLfcCol="log2FoldChange",
                                   direction="both",
                                   chr="all", nameChrCol="chr")

    x<-list(salmon=salmonSig$wormbaseID, Meeuse=oscillating$wormbaseID,
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
            hs_up=hsUP$wormbaseID, hs_down=hsDOWN$wormbaseID)
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



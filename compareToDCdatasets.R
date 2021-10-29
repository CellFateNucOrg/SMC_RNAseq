library(eulerr)
library(lattice)
library(ggplot2)

source("functions.R")
source("./variableSettings.R")
scriptName <- "compareToDCdatasets"
print(scriptName)

if(filterData){
  fileNamePrefix<-filterPrefix
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
} else {
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
}

makeDirs(outPath,dirNameList=paste0(c("plots/"),
                                    paste0(dirname(fileNamePrefix),"/",
                                           scriptName)))


eulerLabelsType<-c("counts")

###########################-
## compare samples -----
##########################-

#######-
## venn diagrams----
#######-

for (grp in useContrasts){
  #grp="dpy26cs"
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

  if(filterData){
    # remove filtered genes
    idx<-salmon$wormbaseID %in% toFilter
    salmon<-salmon[!idx,]
  }

  #####################################################-
  ## compare to public data ----
  #####################################################-

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
    idx<-kramerDpy27$wormbaseID %in% toFilter
    kramerDpy27<-kramerDpy27[!idx,]

    idx<-kramerDpy21$wormbaseID %in% toFilter
    kramerDpy21<-kramerDpy21[!idx,]
  }

  salmonSig<-getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                 namePadjCol="padj",
                                 nameLfcCol="log2FoldChange",
                                 direction="both",
                                 chr="all", nameChrCol="chr")


  DC<-list(salmon=salmonSig$wormbaseID, dpy27=kramerDpy27$wormbaseID,
          dpy21=kramerDpy21$wormbaseID)
  names(DC)<-c(grp, "dpy-27", "dpy-21")
  txtLabels<-list()
  txtLabels[paste0("% ",names(DC)[2]," in ",names(DC)[1])]<-round(100*length(intersect(DC[[1]],DC[[2]]))/length(DC[[2]]),1)
  txtLabels[paste0("% ",names(DC)[3]," in ",names(DC)[1])]<-round(100*length(intersect(DC[[1]],DC[[3]]))/length(DC[[3]]),1)
  fit<-euler(DC)
  percentages<-paste0(txtLabels[[1]], names(txtLabels)[1], " & ",
                      txtLabels[[2]], names(txtLabels)[2])

  pdf(file=paste0(outPath, "/plots/",outputNamePrefix,"venn_",
                  grp,"VsKarmer_padj",
                  padjVal,"_lfc", lfcVal,".pdf"),width=5, height=10, paper="a4")
  p1<-plot(fit, quantities=list(type=eulerLabelsType),
       main=list(label=paste0(grp," vs Kramer(2015): |lfc|>", lfcVal,
                              ", padj<",padjVal,"\n",percentages), fontsize=8, y=0.7))
  print(p1)
  dev.off()

  includeChrX<-ifelse("chrX" %in% salmonSig$chr,T,F)

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

    DC<-list(salmon=salmondc$wormbaseID, dpy27=kramerdpy27dc$wormbaseID,
            dpy21=kramerdpy21dc$wormbaseID)
    names(DC)<-c(grp, "dpy-27", "dpy-21")
    txtLabels<-list()
    txtLabels[paste0("% ",names(DC)[2]," in ",names(DC)[1])]<-round(100*length(intersect(DC[[1]],DC[[2]]))/length(DC[[2]]),1)
    txtLabels[paste0("% ",names(DC)[3]," in ",names(DC)[1])]<-round(100*length(intersect(DC[[1]],DC[[3]]))/length(DC[[3]]),1)

    fit<-euler(DC)
    percentages<-paste0(txtLabels[[1]], names(txtLabels)[1], " & ",
                        txtLabels[[2]], names(txtLabels)[2])

    pdf(file=paste0(outPath, "/plots/",outputNamePrefix,
                    "venn_chrXup_",grp, "VsKramer_padj",
                    padjVal,"_lfc", lfcVal,".pdf"),width=5, height=10, paper="a4")
    p1<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste0("chrX DC genes: lfc>", lfcVal,
                                ", padj<",padjVal,"\n",percentages), fontsize=8, y=0.7))
    print(p1)
    dev.off()
  }

  ###############################-
  ## autosome changed -----
  ###############################-

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

  DC<-list(salmon=salmondc$wormbaseID, dpy27=kramerdpy27dc$wormbaseID,
          dpy21=kramerdpy21dc$wormbaseID)
  names(DC)<-c(grp, "dpy-27", "dpy-21")
  txtLabels<-list()
  txtLabels[paste0("% ",names(DC)[2]," in ",names(DC)[1])]<-round(100*length(intersect(DC[[1]],DC[[2]]))/length(DC[[2]]),1)
  txtLabels[paste0("% ",names(DC)[3]," in ",names(DC)[1])]<-round(100*length(intersect(DC[[1]],DC[[3]]))/length(DC[[3]]),1)

  fit<-euler(DC)
  percentages<-paste0(txtLabels[[1]], names(txtLabels)[1], " & ",
                      txtLabels[[2]], names(txtLabels)[2])

  pdf(file=paste0(outPath, "/plots/",outputNamePrefix,
                  "venn_chrAall_",grp, "VsKramer_padj",
                  padjVal,"_lfc", lfcVal,".pdf"),width=5, height=10,  paper="a4")
  p1<-plot(fit, quantities=list(type=eulerLabelsType),
       main=list(label=paste0("Autosomal genes: |lfc|>", lfcVal,
                              ", padj<",padjVal,"\n",percentages), fontsize=8, y=0.7))
  print(p1)
  dev.off()
}


grp="dpy26cs"
#for(grp in c("dpy26cs")){
for(grp in useContrasts){
  ###############################-
  ## classical dc genes -----
  ###############################-

  pubDC<-readRDS(paste0(outPath,"/publicData/published_DCgr.rds"))
  pubNDC<-readRDS(paste0(outPath,"/publicData/published_NDCgr.rds"))

  if(filterData){
    # remove filtered genes
    idx<-pubDC$wormbaseID %in% toFilter
    pubDC<-pubDC[!idx,]

    idx<-pubNDC$wormbaseID %in% toFilter
    pubNDC<-pubNDC[!idx,]
  }

  # if(includeChrX){
  #   salmondc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
  #   DC<-list(salmon=salmondc$wormbaseID, DC=pubDC$wormbaseID,
  #            nonDC=pubNDC$wormbaseID)
  #   names(DC)<-c(grp, "DC", "nonDC")
  #   txtLabels<-list()
  #   txtLabels[paste0("% ",names(DC)[2]," in ",names(DC)[1])]<-round(100*length(intersect(DC[[1]],DC[[2]]))/length(DC[[2]]),1)
  #   txtLabels[paste0("% ",names(DC)[3]," in ",names(DC)[1])]<-round(100*length(intersect(DC[[1]],DC[[3]]))/length(DC[[3]]),1)
  #
  #   fit<-euler(DC)
  #   percentages<-paste0(txtLabels[[1]], names(txtLabels)[1], " & ",
  #                       txtLabels[[2]], names(txtLabels)[2])
  #
  #   pdf(file=paste0(outPath, "/plots/",outputNamePrefix,"venn_",
  #                   grp, "VsPapers_padj",
  #                   padjVal,"_lfc", lfcVal,".pdf"),width=5, height=10,
  #       paper="a4")
  #   p1<-plot(fit, quantities=list(type=eulerLabelsType),
  #            main=list(label=paste0("chrX dc genes: lfc>", lfcVal, ", padj<",
  #                                   padjVal,"\n",percentages), fontsize=8, y=0.7))
  #   print(p1)
  #   dev.off()
  # }





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

  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

  if(includeChrX){
    salmondc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
    DC<-list(salmon=salmondc$wormbaseID, JansDC=JansDC$wormbaseID,
             JansNDC=JansNDC$wormbaseID)
    names(DC)<-c(grp, "JansDC", "JansNDC")
    txtLabels<-list()
    txtLabels[paste0("% ",names(DC)[2]," in ",names(DC)[1])]<-round(100*length(intersect(DC[[1]],DC[[2]]))/length(DC[[2]]),1)
    txtLabels[paste0("% ",names(DC)[3]," in ",names(DC)[1])]<-round(100*length(intersect(DC[[1]],DC[[3]]))/length(DC[[3]]),1)

    fit<-euler(DC)
    percentages<-paste0(txtLabels[[1]], names(txtLabels)[1], " & ",
                        txtLabels[[2]], names(txtLabels)[2])

    pdf(paste0(outPath, "/plots/",outputNamePrefix,"venn_",grp,
               "VsJans2009_padj", padjVal,"_lfc", lfcVal,".pdf"),
        width=5, height=10, paper="a4")
    p1<-plot(fit, quantities=list(type=eulerLabelsType),
             main=list(label=paste0("chrX dc genes: lfc>", lfcVal, ", padj<",
                                    padjVal,"\n",percentages), fontsize=8,
                       y=0.7))
    print(p1)
    dev.off()
  }
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

  for (grp in useContrasts){
    salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
    salmonSig<-getSignificantGenes(salmon, padj=padjVal,
                                   lfc=lfcVal,
                                   namePadjCol="padj",
                                   nameLfcCol="log2FoldChange",
                                   direction="both",
                                   chr="all", nameChrCol="chr")

    OSC<-list(salmon=salmonSig$wormbaseID, Meeuse=oscillating$wormbaseID,
            Latorre=latorre$wormbaseID)
    names(OSC)<-c(grp, "Meeuse", "Latorre")
    txtLabels<-list()
    txtLabels[paste0("% ",names(OSC)[1]," in ",names(OSC)[2])]<-round(100*length(intersect(OSC[[1]],OSC[[2]]))/length(OSC[[1]]),1)
    txtLabels[paste0("% ",names(OSC)[1]," in ",names(OSC)[3])]<-round(100*length(intersect(OSC[[1]],OSC[[3]]))/length(OSC[[1]]),1)
    txtLabels[paste0("% ",names(OSC)[1]," in union")]<-round(100*length(intersect(OSC[[1]],union(OSC[[2]],OSC[[3]])))/length(OSC[[1]]),1)

    fit<-euler(OSC)
    percentages<-paste0(txtLabels[[1]], names(txtLabels)[1], " &\n",
                        txtLabels[[2]], names(txtLabels)[2], " & ",
                        txtLabels[[3]], names(txtLabels)[3])

    pdf(paste0(outPath, "/plots/",outputNamePrefix,"venn_",grp,"VsOscillating_padj",
               padjVal,"_lfc", lfcVal,".pdf"),width=5, height=10, paper="a4")
    p1<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste0(grp," vs oscillating: |lfc|>", lfcVal,
                                ", padj<",padjVal,"\n",percentages), fontsize=8, y=0.7))
    print(p1)
    dev.off()



    # hsUP
    salmonUp<-getSignificantGenes(salmon, padj=padjVal,
                                   lfc=lfcVal,
                                   namePadjCol="padj",
                                   nameLfcCol="log2FoldChange",
                                   direction="gt",
                                   chr="all", nameChrCol="chr")
    salmonDown<-getSignificantGenes(salmon, padj=padjVal,
                                  lfc= -lfcVal,
                                  namePadjCol="padj",
                                  nameLfcCol="log2FoldChange",
                                  direction="lt",
                                  chr="all", nameChrCol="chr")

    HS<-list(up=salmonUp$wormbaseID, down=salmonDown$wormbaseID,
            hs_up=hsUP$wormbaseID, hs_down=hsDOWN$wormbaseID)
    names(HS)<-c(paste0(grp,c("_up","_down")),
                       "hs_up", "hs_down")
    txtLabels<-list()
    txtLabels[paste0("% ",names(HS)[1]," in ",names(HS)[3])]<-round(100*length(intersect(HS[[1]],HS[[3]]))/length(HS[[1]]),1)
    txtLabels[paste0("% ",names(HS)[1]," in ",names(HS)[4])]<-round(100*length(intersect(HS[[1]],HS[[4]]))/length(HS[[1]]),1)
    txtLabels[paste0("% ",names(HS)[2]," in ",names(HS)[3])]<-round(100*length(intersect(HS[[2]],HS[[3]]))/length(HS[[2]]),1)
    txtLabels[paste0("% ",names(HS)[2]," in ",names(HS)[4])]<-round(100*length(intersect(HS[[2]],HS[[4]]))/length(HS[[2]]),1)

    fit<-euler(HS)
    percentages<-paste0(txtLabels[[1]], names(txtLabels)[1], " & ",
                        txtLabels[[2]], names(txtLabels)[2], " &\n",
                        txtLabels[[3]], names(txtLabels)[3], " & ",
                        txtLabels[[4]], names(txtLabels)[4])
    pdf(paste0(outPath, "/plots/",outputNamePrefix,"venn_",grp,"VsHeatshock_padj",
               padjVal,"_lfc", lfcVal,".pdf"),width=5, height=10, paper="a4")
    p1<-plot(fit, quantities=list(type=eulerLabelsType),
          main=list(label=paste0(grp," vs heatshock: |lfc|>", lfcVal,
                                    ", padj<",padjVal,"\n",percentages), fontsize=8, y=0.7))
    print(p1)
    dev.off()

  }
}



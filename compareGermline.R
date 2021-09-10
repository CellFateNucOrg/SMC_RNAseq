library(ggplot2)
library(EnhancedVolcano)
library(eulerr)
library(lattice)

source("functions.R")
source("./variableSettings.R")
if(filterData){
  fileNamePrefix=filterPrefix
}


eulerLabelsType<-c("counts")

#####################################################-
## compare to germline -soma data-----
#####################################################-

boeck<-read.csv(paste0(outPath,"/publicData/germlineSomaGenes_Boeck2016.csv"),
                stringsAsFactors=F)
reinke<-read.csv(paste0(outPath,"/publicData/germlineGenes_Reinke2004.csv"),
                 stringsAsFactors=F)

if(filterData){
  # remove filtered genes
  idx<-boeck$wormbaseID %in% toFilter
  boeck<-boeck[!idx,]

  idx<-reinke$wormbaseID %in% toFilter
  reinke<-reinke[!idx,]
}


for (grp in useContrasts){
  #grp="dpy26cs"
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

  if(filterData){
    # remove filtered genes
    idx<-salmon$wormbaseID %in% toFilter
    salmon<-salmon[!idx,]
  }

  salmonSig<-getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                 namePadjCol="padj",
                                 nameLfcCol="log2FoldChange",
                                 direction="both",
                                 chr="all", nameChrCol="chr")


  germsoma<-list(salmon=salmonSig$wormbaseID, germline=reinke$wormbaseID)
  names(germsoma)<-c(grp,"germline")

  pdf(file=paste0(outPath, "/plots/",fileNamePrefix,"venn_", grp,
                  "VsGermline_padj", padjVal, "_lfc", lfcVal,".pdf"),
                  width=5,height=10,paper="a4")
  fit<-euler(germsoma)
  totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
              sum(fit$original.values[grep(x, names(fit$original.values))]))}),
                   collapse="  ")
  p1<-plot(fit, quantities=list(type=eulerLabelsType),
           main=list(label=paste0(grp," vs Reinke(2004): |lfc|>", lfcVal,
                                  ", padj<",padjVal,"\n",totalSums), fontsize=8))
  print(p1)

  germsoma<-list(salmon=salmonSig$wormbaseID,
          germlineL4=boeck$wormbaseID[boeck$germline=="germlineL4"],
          somaL4=boeck$wormbaseID[boeck$germline=="somaL4"])
  names(germsoma)<-c(grp,"germlineL4","somaL4")
  txtLabels<-list()
  txtLabels[paste0("% ",names(germsoma)[1]," in ",names(germsoma)[2])]<-round(100*length(intersect(germsoma[[1]],germsoma[[2]]))/length(germsoma[[1]]),1)
  txtLabels[paste0("% ",names(germsoma)[1]," in ",names(germsoma)[3])]<-round(100*length(intersect(germsoma[[1]],germsoma[[3]]))/length(germsoma[[1]]),1)
  fit<-euler(germsoma)
  percentages<-paste0(txtLabels[[1]], names(txtLabels)[1], " & \n",
                      txtLabels[[2]], names(txtLabels)[2])
  p2<-plot(fit, quantities=list(type=eulerLabelsType),
           main=list(label=paste0(grp," vs Boeck(2016): |lfc|>", lfcVal,
                                  ", padj<",padjVal,"\n",percentages), fontsize=8))
  print(p2)


  germsoma<-list(salmon=salmonSig$wormbaseID,
          gonadYA=boeck$wormbaseID[boeck$germline=="gonadYA"],
          somaYA=boeck$wormbaseID[boeck$germline=="somaYA"])
  names(germsoma)<-c(grp,"gonadYA","somaYA")
  txtLabels<-list()
  txtLabels[paste0("% ",names(germsoma)[1]," in ",names(germsoma)[2])]<-round(100*length(intersect(germsoma[[1]],germsoma[[2]]))/length(germsoma[[1]]),1)
  txtLabels[paste0("% ",names(germsoma)[1]," in ",names(germsoma)[3])]<-round(100*length(intersect(germsoma[[1]],germsoma[[3]]))/length(germsoma[[1]]),1)
  fit<-euler(germsoma)
  percentages<-paste0(txtLabels[[1]], names(txtLabels)[1], " & \n",
                      txtLabels[[2]], names(txtLabels)[2])
  p3<-plot(fit, quantities=list(type=eulerLabelsType),
           main=list(label=paste0(grp," vs Boeck(2016): |lfc|>", lfcVal,
                                  ", padj<",padjVal,"\n",percentages), fontsize=8))
  print(p3)
  dev.off()

  #p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
  # ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_",
  #                                 grp,"VsGermline_padj",
  #                                 padjVal,"_lfc", lfcVal,".pdf"),
  #                 plot=p, device="pdf",width=29,height=11,units="cm")
}





## upregulated genes -----
sigTables<-list()
for (grp in useContrasts){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

  sigTables[[grp]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="gt",
                        chr="all", nameChrCol="chr"))
}

sigGenes<-lapply(sigTables, "[[" ,"wormbaseID")
prettySampleNames<-names(sigGenes)
sigGenes[["germline"]]<-reinke$wormbaseID

pdf(file=paste0(outPath, "/plots/",fileNamePrefix,
                "venn_UpVsGermline_padj",
                padjVal,"_lfc", lfcVal,".pdf"),
    width=5,height=10,paper="a4")

eulerPlotList<-list()
for(n in prettySampleNames){
  fit<-euler(sigGenes[c(n,"germline")])
  totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
              sum(fit$original.values[grep(x, names(fit$original.values))]))}),
              collapse="  ")
  eulerPlotList[[n]]<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste0(n, " genes up: lfc>", lfcVal, ", padj<",
                                padjVal,"\n",totalSums), fontsize=8))
}
lapply(eulerPlotList,print)
dev.off()

# p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
# ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,
#                                 "venn_UpVsGermline_padj",
#                                 padjVal,"_lfc", lfcVal,".pdf"),
#                 plot=p, device="pdf",width=29,height=11,units="cm")

if(all(c("kle-2cs","scc-1cs") %in% prettySampleNames)){

  kle2only<-base::setdiff(sigGenes[["kle-2cs"]],sigGenes[["scc-1cs"]])
  scc1only<-base::setdiff(sigGenes[["scc-1cs"]],sigGenes[["kle-2cs"]])

  only<-list(kle2only=kle2only, scc1only=scc1only,germline=sigGenes[["germline"]])
  fit<-euler(only)
  totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
                                                                     sum(fit$original.values[grep(x, names(fit$original.values))]))}),
                   collapse="  ")
  p11<-plot(fit, quantities=list(type=eulerLabelsType),
            main=list(label=paste0("kle-2 or scc-1 genes up: lfc>", lfcVal,
                                   ", padj<",padjVal,"\n",totalSums), fontsize=8))


  kle2scc1<-base::intersect(sigGenes[["kle-2cs"]],sigGenes[["scc-1cs"]])
  both<-list(kle2andscc1=kle2scc1,germline=sigGenes[["germline"]])
  fit<-euler(both)
  totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
                                                                     sum(fit$original.values[grep(x, names(fit$original.values))]))}),
                   collapse="  ")
  p12<-plot(fit, quantities=list(type=eulerLabelsType),
            main=list(label=paste0("kle-2 and scc-1 genes up: lfc>", lfcVal,
                                   ", padj<",padjVal,"\n",totalSums), fontsize=8))
}

## downregulated genes-----
sigTables<-list()
for (grp in useContrasts){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

  sigTables[[grp]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=-lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="lt",
                        chr="all", nameChrCol="chr"))
}

sigGenes<-lapply(sigTables, "[[","wormbaseID")
prettySampleNames<-names(sigGenes)
sigGenes[["germline"]]<-reinke$wormbaseID

pdf(file=paste0(outPath, "/plots/",fileNamePrefix,
                "venn_DownVsGermline_padj",
                padjVal,"_lfc", lfcVal,".pdf"),
    width=5,height=10,paper="a4")



eulerPlotList<-list()
for(n in prettySampleNames){
  fit<-euler(sigGenes[c(n,"germline")])
  totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
                                                                     sum(fit$original.values[grep(x, names(fit$original.values))]))}),
                   collapse="  ")
  eulerPlotList[[n]]<-plot(fit, quantities=list(type=eulerLabelsType),
           main=list(label=paste0(n, " genes down: lfc < -", lfcVal,
                                  ", padj<",padjVal,"\n",totalSums), fontsize=8))
}
lapply(eulerPlotList,print)
dev.off()

#p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
#ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,
#                                "venn_DownVsGermline_padj",
#                                padjVal,"_lfc", lfcVal,".pdf"),
#                plot=p, device="pdf",width=29,height=11,units="cm")


if(all(c("kle-2cs","scc-1cs") %in% prettySampleNames)){
  kle2only<-base::setdiff(sigGenes[["kle-2cs"]],sigGenes[["scc-1cs"]])
  scc1only<-base::setdiff(sigGenes[["scc-1cs"]],sigGenes[["kle-2cs"]])

  only<-list(kle2only=kle2only, scc1only=scc1only,germline=sigGenes[["germline"]])
  fit<-euler(only)
  totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
                                                                     sum(fit$original.values[grep(x, names(fit$original.values))]))}),
                   collapse="  ")
  p13<-plot(fit, quantities=list(type=eulerLabelsType),
            main=list(label=paste0("kle-2 or scc-1 genes down: lfc>", lfcVal,
                                   ", padj<",padjVal,"\n",totalSums), fontsize=8))

  kle2scc1<-base::intersect(sigGenes[["kle-2cs"]],sigGenes[["scc-1cs"]])

  both<-list(kle2andscc1=kle2scc1,germline=sigGenes[["germline"]])
  fit<-euler(both)
  totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
                                                                     sum(fit$original.values[grep(x, names(fit$original.values))]))}),
                   collapse="  ")
  p14<-plot(fit, quantities=list(type=eulerLabelsType),
            main=list(label=paste0("kle-2 and scc-1 genes down: lfc>", lfcVal,
                                   ", padj<",padjVal,"\n",totalSums), fontsize=8))

  pdf(file=paste0(outPath, "/plots/",fileNamePrefix,
                  "venn_kle2scc1setsVsGermline_padj",
                  padjVal,"_lfc", lfcVal,".pdf"),
      width=5,height=10,paper="a4")
  print(p11)
  print(p12)
  print(p13)
  print(p14)
  dev.off()
}
# p<-ggpubr::ggarrange(p11,p12,p13,p14,ncol=2,nrow=2)
# ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,
#                                 "venn_kle2scc1setsVsGermline_padj",
#                                 padjVal,"_lfc", lfcVal,".pdf"),
#                 plot=p, device="pdf",width=21,height=19,units="cm")

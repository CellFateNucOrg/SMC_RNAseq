library(ggplot2)
library(EnhancedVolcano)
library(eulerr)
library(lattice)
library(fgsea)
library(clusterProfiler)

source("functions.R")
source("./variableSettings.R")
scriptName <- "compareGermline"
print(scriptName)

if(filterData){
  fileNamePrefix<-filterPrefix
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
} else {
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
}

makeDirs(outPath,dirNameList=paste0(c("plots/","txt/"),
                                    paste0(dirname(fileNamePrefix),"/",
                                           scriptName)))

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

  pdf(file=paste0(outPath, "/plots/",outputNamePrefix,"venn_", grp,
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
  # ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,"venn_",
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

pdf(file=paste0(outPath, "/plots/",outputNamePrefix,
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
# ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
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

pdf(file=paste0(outPath, "/plots/",outputNamePrefix,
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
#ggplot2::ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,
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

  pdf(file=paste0(outPath, "/plots/",outputNamePrefix,
                  "venn_kle2scc1setsVsGermline_padj",
                  padjVal,"_lfc", lfcVal,".pdf"),
      width=5,height=10,paper="a4")
  print(p11)
  print(p12)
  print(p13)
  print(p14)
  dev.off()
}



#####################################################-
## GSEA germline-soma data-----
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

boeckLst<-split(boeck$wormbaseID,boeck$germline)
names(boeckLst)<-paste0("Boeck2016_",names(boeckLst))

reinkeLst<-split(reinke$wormbaseID,reinke$exclusive.category)
names(reinkeLst)<-paste0("Reinke2004_",gsub(" ",".",names(reinkeLst)))

germline<-c(boeckLst,reinkeLst)

gseaTbl<-list()
leadEdgeTbl<-NULL
for (grp in useContrasts){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

  if(filterData){
    # remove filtered genes
    idx<-salmon$wormbaseID %in% toFilter
    salmon<-salmon[!idx,]
  }
  # making ranks (ranks go from high to low)
  ranks <- salmon$log2FoldChange
  names(ranks) <- salmon$wormbaseID
  head(ranks)

  fgseaRes <- fgsea(germline, ranks, minSize=5, maxSize = 5000)

  head(fgseaRes[order(padj), ])
  fgseaRes[padj<=0.05,]
  gseaTbl[[grp]]<-plotGseaTable(germline, ranks, fgseaRes,gseaParam=0.1,render=F,
                                colwidths = c(5, 3, 0.8, 0, 1.2))
  #barplot(sort(ranks, decreasing = T))
  gseaList<-list()
  for (pathName in unlist(fgseaRes[padj<0.05,"pathway"])) {
    gseaList[[pathName]]<-plotEnrichment(germline[[pathName]], ranks) +
      labs(title=paste0(grp," enrichment in ",pathName)) +
      geom_vline(xintercept=sum(sort(ranks)>0), colour="grey40")+
      annotate("text", x=length(ranks)*0.75, y=fgseaRes[pathway==pathName,ES]*0.9,
               label= paste(paste(paste(c("padj","NES","size"),
                                        fgseaRes[pathway==pathName,c(round(padj,3),round(NES,1),size)],
                                        sep=":"),collapse=", ")))
    leadEdge<-unlist(fgseaRes[fgseaRes$pathway==pathName,"leadingEdge"])
    leadEdge<-data.frame(wormbaseID=leadEdge)
    #row.names(leadEdge)<-NULL
    leadEdge<-left_join(leadEdge,as.data.frame(salmon), by="wormbaseID")
    leadEdge$group<-grp
    leadEdge$pathway<-pathName
    if(is.null(leadEdgeTbl)){
      leadEdgeTbl<-leadEdge
    } else {
      leadEdgeTbl<-rbind(leadEdgeTbl,leadEdge)
    }
    #print(paste0(grp," in ",pathName,":  ",paste(sort(leadEdge$publicID),collapse=",")))
    #print(leadEdge)
  }
  pdf(file=paste0(outPath, "/plots/",outputNamePrefix,"gsea_", grp,
                  "VsGermlineDatasets.pdf"),
      width=5,height=10,paper="a4")
  #p<-ggpubr::ggarrange(plotlist=gseaList,nrow=2,ncol=1)
  p<-gridExtra::marrangeGrob(grobs=gseaList, nrow=3, ncol=1)
  print(p)
  dev.off()
}

if(length(gseaList)>0){
  pdf(file=paste0(outPath, "/plots/",outputNamePrefix,"gseaAll_GermlineDatasets.pdf"),
      width=16,height=11)
  p<-gridExtra::marrangeGrob(grobs=gseaTbl, ncol=3, nrow=1,padding=unit(0.01,"line"))
  print(p)
  dev.off()

  write.table(leadEdgeTbl,file=paste0(outPath,"/txt/",outputNamePrefix,"gseaLeadEdgeGenes.tsv"),sep="\t")
}

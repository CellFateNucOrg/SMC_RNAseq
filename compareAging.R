library(ggplot2)
library(eulerr)
library(lattice)
library(dplyr)
library(tibble)
library(fgsea)
library(clusterProfiler)
#library(plotenrich)

source("functions.R")
source("./variableSettings.R")
scriptName <- "compareAging"
print(scriptName)

if(filterData){
  fileNamePrefix<-filterPrefix
  # #add hs gene filteration
  # hsUp<-readRDS(paste0(outPath,"/publicData/hsUp_garrigues2019.rds"))
  # hsDown<-readRDS(paste0(outPath,"/publicData/hsDown_garrigues2019.rds"))
  # toFilter<-unique(sort(c(toFilter,hsUp$wormbaseID,hsDown$wormbaseID)))
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
} else {
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
}

makeDirs(outPath,dirNameList=paste0(c("plots/","txt/"),
                                    paste0(dirname(fileNamePrefix),"/",
                                           scriptName)))

#eulerLabelsType<-c("counts")


# #####################################################-
# ## compare to age regulated genes (Budovskaya 2008)-----
# #####################################################-
#
# ageReg<-read.csv(paste0(outPath,"/publicData/AgeRegulated_Budovskaya2008.csv"),
#                 stringsAsFactors=F)
#
# tcMA<-read.csv(paste0(outPath,"/publicData/AgingTC_Budovskaya2008.csv"),
#                  stringsAsFactors=F)
#
# localPval=0.05
# localLFCval=0
# age1MA<-read.delim(paste0(outPath,"/publicData/age1MA_Budovskaya2008.tsv"),
#                  stringsAsFactors=F)
# age1MA$description<-NULL
#
# daf16MA<-read.delim(paste0(outPath,"/publicData/daf16MA_Budovskaya2008.tsv"),
#                    stringsAsFactors=F)
# daf16MA$description<-NULL
#
# if(filterData){
#   # remove filtered genes
#   idx<-ageReg$wormbaseID %in% toFilter
#   ageReg<-ageReg[!idx,] #1011 down to 583 genes
#   idx<-tcMA$wormbaseID %in% toFilter
#   tcMA<-tcMA[!idx,] #15752 down to 12270 genes
#   idx<-age1MA$wormbaseID %in% toFilter
#   age1MA<-age1MA[!idx,] #14645 down to 11400 genes
#   idx<-daf16MA$wormbaseID %in% toFilter
#   daf16MA<-daf16MA[!idx,] #14645 down to 11400 genes
#
# }
#
# ageRegUp<-ageReg[ageReg$day11_avg>0,]
# ageRegDown<-ageReg[ageReg$day11_avg<0,]
#
# age1MAup<-age1MA[age1MA$pval<localPval & age1MA$ratio > localLFCval,] #851
# daf16MAdown<-daf16MA[daf16MA$pval<localPval & daf16MA$ratio < localLFCval,] #2909
#
# for (grp in useContrasts){
#   #grp="dpy26cs"
#   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
#
#   if(filterData){
#     # remove filtered genes
#     idx<-salmon$wormbaseID %in% toFilter
#     salmon<-salmon[!idx,]
#   }
#
#   # salmonSig<-getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
#   #                                namePadjCol="padj",
#   #                                nameLfcCol="log2FoldChange",
#   #                                direction="both",
#   #                                chr="all", nameChrCol="chr")
#
#   salmonSigUp<-getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
#                                    namePadjCol="padj",
#                                    nameLfcCol="log2FoldChange",
#                                    direction="gt",
#                                    chr="all", nameChrCol="chr")
#
#   salmonSigDown<-getSignificantGenes(salmon, padj=padjVal, lfc=-lfcVal,
#                                      namePadjCol="padj",
#                                      nameLfcCol="log2FoldChange",
#                                      direction="lt",
#                                      chr="all", nameChrCol="chr")
#
#
#   pdf(file=paste0(outPath, "/plots/",outputNamePrefix,"venn_", grp,
#                   "VsBudovskaya2008_padj", padjVal, "_lfc", lfcVal,".pdf"),
#       width=5,height=10,paper="a4")
#
#   geneList<-list(salmon=salmonSigUp$wormbaseID, ageRegUp=ageRegUp$wormbaseID,
#                  ageRegDown=ageRegDown$wormbaseID)
#   names(geneList)<-c(paste0(grp," up"),"Age upregulated", "Age downregulated")
#
#   fit<-euler(geneList)
#   totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
#                                                                      sum(fit$original.values[grep(x, names(fit$original.values))]))}),
#                    collapse="  ")
#   p1<-plot(fit, quantities=list(type=eulerLabelsType),
#            main=list(label=paste0(grp,"up vs Budovskaya(2008): |lfc|>", lfcVal,
#                                   ", padj<",padjVal,"\n",totalSums), fontsize=8))
#   print(p1)
#
#   geneList<-list(salmon=salmonSigDown$wormbaseID, ageRegUp=ageRegUp$wormbaseID,
#                  ageRegDown=ageRegDown$wormbaseID)
#   names(geneList)<-c(paste0(grp," down"),"day11 up", "day11 down")
#
#   fit<-euler(geneList)
#   totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
#               sum(fit$original.values[grep(x, names(fit$original.values))]))}),
#                    collapse="  ")
#   p1a<-plot(fit, quantities=list(type=eulerLabelsType),
#            main=list(label=paste0(grp," down vs Budovskaya(2008): |lfc|>", lfcVal,
#                                   ", padj<",padjVal,"\n",totalSums), fontsize=8))
#   print(p1a)
#
#   # compare to age-1 up and daf-16 down genes
#   geneList<-list(salmon=salmonSigUp$wormbaseID,age1=age1MAup$wormbaseID ,daf16=daf16MAdown$wormbaseID)
#   names(geneList)<-c(paste0(grp," up"),"age-1 up","daf-16 down")
#
#   fit<-euler(geneList)
#   totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
#               sum(fit$original.values[grep(x, names(fit$original.values))]))}),
#                    collapse="  ")
#   p2<-plot(fit, quantities=list(type=eulerLabelsType),
#            main=list(label=paste0(grp," up vs Budovskaya(2008): |lfc|>", lfcVal,
#                                   ", padj<",padjVal,"\n",totalSums), fontsize=8))
#   print(p2)
#
#   geneList<-list(salmon=salmonSigDown$wormbaseID,age1=age1MAup$wormbaseID ,daf16=daf16MAdown$wormbaseID)
#   names(geneList)<-c(paste0(grp," down"),"age-1 up","daf-16 down")
#
#   fit<-euler(geneList)
#   totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
#               sum(fit$original.values[grep(x, names(fit$original.values))]))}),
#                                       collapse="  ")
#   p3<-plot(fit, quantities=list(type=eulerLabelsType),
#                               main=list(label=paste0(grp," down vs Budovskaya(2008): |lfc|>", lfcVal,
#                                                      ", padj<",padjVal,"\n",totalSums), fontsize=8))
#   print(p3)
#
#   idx<-salmonSigUp$wormbaseID %in% ageRegDown$wormbaseID |
#            salmonSigUp$wormbaseID %in% age1MAup$wormbaseID |
#            salmonSigUp$wormbaseID %in% daf16MAdown$wormbaseID
#   shared<-as.data.frame(salmonSigUp[idx,])
#   shared<-left_join(shared,data.frame(tcMA[,c("wormbaseID","day11_avg")]),by="wormbaseID")
#   shared<-left_join(shared,data.frame(wormbaseID=age1MA$wormbaseID,age1MA=age1MA$ratio),by="wormbaseID")
#   shared<-left_join(shared,data.frame(wormbaseID=daf16MA$wormbaseID,daf16MA=daf16MA$ratio),by="wormbaseID")
#   write.table(shared,file=paste0(outPath,"/txt/",outputNamePrefix,"aging_Budovskaya2008_",grp,
#                                  "-up_padj", padjVal, "_lfc", lfcVal,".tsv"),sep="\t")
#   par(mfrow=c(4,1))
#
#   lowerFn <- function(data, mapping, method = "lm", ...) {
#     p <- ggplot(data = data, mapping = mapping) +
#       geom_point(colour = "blue") +
#       geom_smooth(method = method, color = "red", ...)
#     p
#   }
#   p4<-GGally::ggpairs(data=shared,columns=c("log2FoldChange","day11_avg","age1MA","daf16MA"),
#                       lower=list(continuous = GGally::wrap(lowerFn, method = "lm"))) +
#     ggtitle(paste0(grp))
#   print(p4)
#   dev.off()
# }
#
#
# #####################################################-
# ## compare to age regulated genes (Murphy 2003)-----
# #####################################################-
#
# ageClassI<-read.csv(paste0(outPath,"/publicData/agingClassI_Murphy2003.csv"),
#                  stringsAsFactors=F,sep=";")
# ageClassI$description<-NULL
# ageClassII<-read.csv(paste0(outPath,"/publicData/agingClassII_Murphy2003.csv"),
#                     stringsAsFactors=F,sep=";")
# ageClassII$description<-NULL
#
# if(filterData){
#   # remove filtered genes
#   idx<-ageClassI$wormbaseID %in% toFilter
#   ageClassI<-ageClassI[!idx,] #237 down to 132 genes
#   idx<-ageClassII$wormbaseID %in% toFilter
#   ageClassII<-ageClassII[!idx,] #231 down to 159 genes
# }
#
# for (grp in useContrasts){
#   #grp="dpy26cs"
#   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
#
#   if(filterData){
#     # remove filtered genes
#     idx<-salmon$wormbaseID %in% toFilter
#     salmon<-salmon[!idx,]
#   }
#
#   salmonSigUp<-getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
#                                    namePadjCol="padj",
#                                    nameLfcCol="log2FoldChange",
#                                    direction="gt",
#                                    chr="all", nameChrCol="chr")
#
#   salmonSigDown<-getSignificantGenes(salmon, padj=padjVal, lfc=-lfcVal,
#                                      namePadjCol="padj",
#                                      nameLfcCol="log2FoldChange",
#                                      direction="lt",
#                                      chr="all", nameChrCol="chr")
#
#
#   pdf(file=paste0(outPath, "/plots/",outputNamePrefix,"venn_", grp,
#                   "VsMurphy2003_padj", padjVal, "_lfc", lfcVal,".pdf"),
#       width=5,height=10,paper="a4")
#
#   geneList<-list(salmon=salmonSigUp$wormbaseID, agingClassI=ageClassI$wormbaseID,
#                  agingClassII=ageClassII$wormbaseID)
#   names(geneList)<-c(paste0(grp," up"),"Aging ClassI", "Aging ClassII")
#
#   fit<-euler(geneList)
#   totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
#                                                                      sum(fit$original.values[grep(x, names(fit$original.values))]))}),
#                    collapse="  ")
#   p1<-plot(fit, quantities=list(type=eulerLabelsType),
#            main=list(label=paste0(grp,"up vs Murphy(2003): |lfc|>", lfcVal,
#                                   ", padj<",padjVal,"\n",totalSums), fontsize=8))
#   print(p1)
#
#   geneList<-list(salmon=salmonSigDown$wormbaseID, agingClassI=ageClassI$wormbaseID,
#                  agingClassII=ageClassII$wormbaseID)
#   names(geneList)<-c(paste0(grp," down"),"Aging ClassI", "Aging ClassII")
#
#   fit<-euler(geneList)
#   totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
#                                                                      sum(fit$original.values[grep(x, names(fit$original.values))]))}),
#                    collapse="  ")
#   p1a<-plot(fit, quantities=list(type=eulerLabelsType),
#             main=list(label=paste0(grp," down vs Murphy(2003): |lfc|>", lfcVal,
#                                    ", padj<",padjVal,"\n",totalSums), fontsize=8))
#   print(p1a)
#
#   salmonSigUp$class<-NA
#   idx<-salmonSigUp$wormbaseID %in% ageClassI$wormbaseID
#   salmonSigUp$class[idx]<-"classI"
#   idx<-salmonSigUp$wormbaseID %in% ageClassII$wormbaseID
#   salmonSigUp$class[idx]<-"classII"
#   shared<-as.data.frame(salmonSigUp[!is.na(salmonSigUp$class),])
#   write.table(shared,file=paste0(outPath,"/txt/",outputNamePrefix,"aging_Murphy2003_",grp,
#                                  "-up_padj", padjVal, "_lfc", lfcVal,".tsv"),sep="\t")
#
#   dev.off()
# }
#
#
#
# #####################################################-
# ## compare to age regulated genes Tarkhov(2019)-----
# #####################################################-
#
# downSignatr<-read.csv(paste0(outPath,"/publicData/agingDownSignature_Tarkhov2019.csv"),
#                     stringsAsFactors=F,sep=";")
# downSignatr$description<-NULL
# upSignatr<-read.csv(paste0(outPath,"/publicData/agingUpSignature_Tarkhov2019.csv"),
#                      stringsAsFactors=F,sep=";")
# upSignatr$description<-NULL
#
# ageBiomrkr<-read.csv(paste0(outPath,"/publicData/agingBiomarkers_Tarkhov2019.csv"),
#                      stringsAsFactors=F,sep=";")
# ageBiomrkr$description<-NULL
#
# if(filterData){
#   # remove filtered genes
#   idx<-downSignatr$wormbaseID %in% toFilter
#   downSignatr<-downSignatr[!idx,] #67 down to 32 genes
#   idx<-upSignatr$wormbaseID %in% toFilter
#   upSignatr<-upSignatr[!idx,] #260 down to 169
#   idx<-ageBiomrkr$wormbaseID %in% toFilter
#   ageBiomrkr<-ageBiomrkr[!idx,] #71 down to 38 genes
# }
#
# upBiomrkr<-ageBiomrkr[ageBiomrkr$Coefficient>0,]
# downBiomrkr<-ageBiomrkr[ageBiomrkr$Coefficient<0,]
#
# for (grp in useContrasts){
#   #grp="dpy26cs"
#   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
#
#   if(filterData){
#     # remove filtered genes
#     idx<-salmon$wormbaseID %in% toFilter
#     salmon<-salmon[!idx,]
#   }
#
#   salmonSigUp<-getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
#                                    namePadjCol="padj",
#                                    nameLfcCol="log2FoldChange",
#                                    direction="gt",
#                                    chr="all", nameChrCol="chr")
#
#   salmonSigDown<-getSignificantGenes(salmon, padj=padjVal, lfc=-lfcVal,
#                                      namePadjCol="padj",
#                                      nameLfcCol="log2FoldChange",
#                                      direction="lt",
#                                      chr="all", nameChrCol="chr")
#
#
#   pdf(file=paste0(outPath, "/plots/",outputNamePrefix,"venn_", grp,
#                   "VsTarkhov2019_padj", padjVal, "_lfc", lfcVal,".pdf"),
#       width=5,height=10,paper="a4")
#
#   geneList<-list(salmon=salmonSigUp$wormbaseID, upSignature=upSignatr$wormbaseID,
#                  downSignature=downSignatr$wormbaseID)
#   names(geneList)<-c(paste0(grp," up"),"Aging up signature", "Aging down signature")
#
#   fit<-euler(geneList)
#   totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
#                                                                      sum(fit$original.values[grep(x, names(fit$original.values))]))}),
#                    collapse="  ")
#   p1<-plot(fit, quantities=list(type=eulerLabelsType),
#            main=list(label=paste0(grp,"up vs Tarkhov(2019): |lfc|>", lfcVal,
#                                   ", padj<",padjVal,"\n",totalSums), fontsize=8))
#   print(p1)
#
#   geneList<-list(salmon=salmonSigDown$wormbaseID, upSignature=upSignatr$wormbaseID,
#                  downSignature=downSignatr$wormbaseID)
#   names(geneList)<-c(paste0(grp," down"),"Aging up signature", "Aging down signature")
#
#   fit<-euler(geneList)
#   totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
#                                                                      sum(fit$original.values[grep(x, names(fit$original.values))]))}),
#                    collapse="  ")
#   p1a<-plot(fit, quantities=list(type=eulerLabelsType),
#             main=list(label=paste0(grp," down vs Tarkhov(2019): |lfc|>", lfcVal,
#                                    ", padj<",padjVal,"\n",totalSums), fontsize=8))
#   print(p1a)
#
#
#   salmonSigUp$class<-NA
#   idx<-salmonSigUp$wormbaseID %in% upSignatr$wormbaseID
#   salmonSigUp$class[idx]<-"upSignatr"
#   idx<-salmonSigUp$wormbaseID %in% downSignatr$wormbaseID
#   salmonSigUp$class[idx]<-"downSignatr"
#   shared<-as.data.frame(salmonSigUp[!is.na(salmonSigUp$class),])
#   write.table(shared,file=paste0(outPath,"/txt/",outputNamePrefix,"aging_Tarkhov2019_",grp,
#                                  "-up_padj", padjVal, "_lfc", lfcVal,".tsv"),sep="\t")
#
#   lowerFn <- function(data, mapping, method = "lm", ...) {
#     p <- ggplot(data = data, mapping = mapping) +
#       geom_point(colour = "blue") +
#       geom_smooth(method = method, color = "red", ...)
#     p
#   }
#
#   forCor<-inner_join(ageBiomrkr,as.data.frame(salmonSig),by="wormbaseID")
#   # only three points
#
#   dev.off()
# }






#####################################################-
## Gene set enrichment - aging time course -----
#####################################################-
# note: gseaplot2 function from enrichplot function supports multiple gene sets
# in a single plot:
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html?q=enrichplot#bar-plot

gseaTbl<-list()
leadEdgeTbl<-NULL
for (grp in useContrasts){
  #grp="dpy26"
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

  agingTCdata<-readRDS(paste0(outPath,"/publicData/agingTCdata.RDS"))

  fgseaRes <- fgsea(agingTCdata, ranks, minSize=5, maxSize = 5000)

  fgseaRes[order(padj), ]
  fgseaRes[padj<=0.05,]
  gseaTbl[[grp]]<-plotGseaTable(agingTCdata, ranks, fgseaRes,gseaParam=0.1,render=F,
                                colwidths = c(5, 3, 0.8, 0, 1.2))
  #barplot(sort(ranks, decreasing = T))
  gseaList<-list()
  for (pathName in unlist(fgseaRes[padj<0.05,"pathway"])) {
    gseaList[[pathName]]<-plotEnrichment(agingTCdata[[pathName]], ranks) +
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
                  "VsAgingTCdata.pdf"),
      width=5,height=10,paper="a4")
  #p<-ggpubr::ggarrange(plotlist=gseaList,nrow=2,ncol=1)
  p<-gridExtra::marrangeGrob(grobs=gseaList, nrow=3, ncol=1)
  print(p)
  dev.off()
}

pdf(file=paste0(outPath, "/plots/",outputNamePrefix,"gseaAll_AgingTCdata.pdf"),
    width=16,height=11)
p<-gridExtra::marrangeGrob(grobs=gseaTbl, ncol=3, nrow=1,padding=unit(0.01,"line"))
print(p)
dev.off()


write.table(leadEdgeTbl,file=paste0(outPath,"/txt/",outputNamePrefix,"gseaLeadEdge_agingTCdata.tsv"),sep="\t",row.names=F,
            quote=F)


#####################################################-
## Gene set enrichment - aging regulators -----
#####################################################-
# note: gseaplot2 function from enrichplot function supports multiple gene sets
# in a single plot:
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html?q=enrichplot#bar-plot

gseaTbl<-list()
leadEdgeTbl<-NULL
for (grp in useContrasts){
  #grp="dpy26"
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

  agingRegData<-readRDS(paste0(outPath,"/publicData/agingRegData.RDS"))

  fgseaRes <- fgsea(agingRegData, ranks, minSize=5, maxSize = 5000)

  fgseaRes[order(padj), ]
  fgseaRes[padj<=0.05,]
  gseaTbl[[grp]]<-plotGseaTable(agingRegData, ranks, fgseaRes,gseaParam=0.1,render=F,
                                colwidths = c(5, 3, 0.8, 0, 1.2))
  #barplot(sort(ranks, decreasing = T))
  gseaList<-list()
  for (pathName in unlist(fgseaRes[padj<0.05,"pathway"])) {
    gseaList[[pathName]]<-plotEnrichment(agingRegData[[pathName]], ranks) +
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
                  "VsAgingRegData.pdf"),
      width=5,height=10,paper="a4")
  #p<-ggpubr::ggarrange(plotlist=gseaList,nrow=2,ncol=1)
  p<-gridExtra::marrangeGrob(grobs=gseaList, nrow=3, ncol=1)
  print(p)
  dev.off()
}

pdf(file=paste0(outPath, "/plots/",outputNamePrefix,"gseaAll_AgingRegData.pdf"),
    width=16,height=11)
p<-gridExtra::marrangeGrob(grobs=gseaTbl, ncol=3, nrow=1,padding=unit(0.01,"line"))
print(p)
dev.off()


write.table(leadEdgeTbl,file=paste0(outPath,"/txt/",outputNamePrefix,"gseaLeadEdge_agingRegData.tsv"),sep="\t",row.names=F,
            quote=F)


keepPath<-c("daf16down-daf2BG_Riedel2013","daf2up_Riedel2013","daf2up_Zarse2012")
leadEdgeTbl<-leadEdgeTbl[leadEdgeTbl$pathway %in% keepPath,]

keepGroup<-c("aux_sdc3BG","dpy26")
leadEdgeTbl<-leadEdgeTbl[leadEdgeTbl$group %in% keepGroup,]
dim(leadEdgeTbl) #1573

leadEdgeTbl %>% dplyr::count(wormbaseID) #895 unique genes
leadEdgeTbl %>% dplyr::group_by(group) %>% dplyr::count(wormbaseID) #1142
leadEdgeTbl %>% dplyr::group_by(pathway) %>% dplyr::count(wormbaseID) #1243
leadEdgeTbl %>% dplyr::count(wormbaseID) %>% dplyr::filter(n>1) #456 genes appear more than once

# in at least two paths or 2 groups
df<-leadEdgeTbl %>% dplyr::group_by(wormbaseID) %>% dplyr::summarise(numPaths=n_distinct(pathway),numGroups=n_distinct(group)) %>% dplyr::filter(numGroups==2)

df1<-left_join(df,leadEdgeTbl[,c("wormbaseID","publicID","sequenceID","chr")],
               by="wormbaseID") %>% distinct() %>% arrange(chr) %>% print(n=Inf)
#247
df1[df1$chr=="chrX",]#166

# in at least two paths or 2 groups - 433 genes
df<-leadEdgeTbl %>% dplyr::group_by(wormbaseID) %>% dplyr::summarise(numPaths=n_distinct(pathway),numGroups=n_distinct(group)) %>% dplyr::filter(numGroups==2 & numPaths==3)

df1<-left_join(df,leadEdgeTbl[,c("wormbaseID","publicID","sequenceID","chr")],
               by="wormbaseID") %>% distinct() %>% arrange(chr) %>% print(n=Inf) #23
df1[df1$chr=="chrX",]#9

write.table(df1,paste0(outPath,"/txt/",outputNamePrefix,"gseaLeadEdge_allPath&Group_agingRegData.tsv"),sep="\t",row.names=F,
            quote=F)

# only on X
df2<-leadEdgeTbl %>% dplyr::filter(chr=="chrX") %>% dplyr::group_by(wormbaseID) %>% dplyr::summarise(numPaths=n_distinct(pathway),numGroups=n_distinct(group))%>% dplyr::filter(numGroups==2)

df3<-left_join(df2,leadEdgeTbl[,c("wormbaseID","publicID","sequenceID","chr")],
               by="wormbaseID") %>% distinct() %>% print(n=Inf) # 166 genes

write.table(df3,paste0(outPath,"/txt/",outputNamePrefix,"gseaLeadEdge_chrXallGroup_agingRegData.tsv"),sep="\t",row.names=F,
            quote=F)



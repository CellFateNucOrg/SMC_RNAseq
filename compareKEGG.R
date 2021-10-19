library(dplyr)
library(tibble)
library(clusterProfiler)
library(pathview)
library(enrichplot)

source("functions.R")
source("./variableSettings.R")

if(filterData){
  fileNamePrefix<-filterPrefix
}

#########################-
# KEGG pathway enrichment -----
#########################-

makeDirs(outPath,dirNameList=paste0("kegg/p",padjVal,"_lfc",lfcVal))

#search_kegg_organism('cel', by='kegg_code')

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

  ###### TODO: need to fix code to be automated
  ranks<-salmonSig$log2FoldChange
  names(ranks)<-paste0("CELE_",salmonSig$sequenceID)
  ranks<-sort(ranks,decreasing=T)

  kk <- enrichKEGG(gene         = unique(names(ranks)),
                   organism     = 'cel',
                   pvalueCutoff = 0.05)
  #head(kk)[c(1:7,9)]
  if(dim(kk)[1]>0){
    write.table(kk,file=paste0(outPath,"/kegg/",fileNamePrefix,"enrichKegg_",grp,
                                  "_padj", padjVal, "_lfc", lfcVal,".tsv"),
                sep="\t",row.names=F)
    # get global plots of pathways
    plotList<-list()
    plotList[["dotplot"]]<-dotplot(kk, showCategory=20) + ggtitle(grp)
    plotList[["cnetstar"]] <- cnetplot(kk, foldChange=ranks, colorEdge = TRUE,
                                       showCategory=10)
    plotList[["cnetcircle"]] <- cnetplot(kk, foldChange=ranks, circular = TRUE,
                                         colorEdge = TRUE, showCategory=10)
    plotList[["heatplot"]] <- heatplot(kk, foldChange=ranks, showCategory=30)

    k2 <- pairwise_termsim(kk)
    if(dim(k2)[1]>6){
      plotList[["treeplot"]] <- treeplot(k2, fontsize=5)
    }
    plotList[["emapplot"]] <- emapplot(k2, cex_category=1,
                                       cex_label_category=0.6, cex_line=0.5) +
      ggtitle(grp)
    plotList[["upsetplot"]] <- upsetplot(kk, n=10)

    pdf(file=paste0(outPath, "/kegg/",fileNamePrefix,"kegg_ORA_", grp,".pdf"),
        width=11,height=8,paper="a4r")
    #p<-ggpubr::ggarrange(plotlist=gseaList,nrow=2,ncol=1)
    p<-gridExtra::marrangeGrob(grobs=plotList, nrow=1, ncol=1)
    print(p)
    dev.off()
    # get individual kegg pathways
    for(pathID in kk$ID){
      if(pathID!="cel01100"){
        pathview(gene.data=ranks, pathway.id=pathID, gene.idtype="kegg",
               kegg.dir=paste0(outPath,"/kegg/p",padjVal,"_lfc",lfcVal),
               out.suffix=paste0("ORA_",grp),
               species="cel")
        system(paste0("mv ",outPath,"/",pathID,".ORA_",grp,".png ",
                      outPath,"/kegg/p",padjVal,"_lfc",lfcVal,"/"))
      }
    }
  }

  kk2 <- gseKEGG(geneList     = ranks[!duplicated(names(ranks))],
               organism     = 'cel',
               minGSSize    = 20,
               pvalueCutoff = 0.05,
               verbose      = TRUE)
  #head(kk2)[c(1:10)]

  if(dim(kk2)[1]>0){
    write.table(kk2,file=paste0(outPath,"/kegg/",fileNamePrefix,"gseKegg_",grp,
                               "_padj", padjVal, "_lfc", lfcVal,".tsv"),
                sep="\t",row.names=F)
    plotList<-list()
    plotList[["dotplot"]]<-dotplot(kk2, showCategory=20) + ggtitle(grp)
    plotList[["cnetstar"]] <- cnetplot(kk2, foldChange=ranks, colorEdge = TRUE,
                                       showCategory=10)
    plotList[["cnetcircle"]] <- cnetplot(kk2, foldChange=ranks, circular = TRUE,
                                         colorEdge = TRUE, showCategory=10)
    plotList[["heatplot"]] <- heatplot(kk2, foldChange=ranks, showCategory=30)

    k22 <- pairwise_termsim(kk2)
    if(dim(k22)[1]>6){
      plotList[["treeplot"]] <- treeplot(k22, hclust_method="average", fontsize=5)
    }
    plotList[["emapplot"]] <- emapplot(k22, cex_category=1,
                                       cex_label_category=0.6, cex_line=0.5) +
      ggtitle(grp)
    plotList[["upsetplot"]] <- upsetplot(kk2, n=10)

    pdf(file=paste0(outPath, "/kegg/",fileNamePrefix,"kegg_GSEA_", grp,".pdf"),
        width=11,height=8,paper="a4r")
    #p<-ggpubr::ggarrange(plotlist=gseaList,nrow=2,ncol=1)
    p<-gridExtra::marrangeGrob(grobs=plotList, nrow=1, ncol=1)
    print(p)
    dev.off()
    for(pathID in kk2$ID){
      if(pathID!="cel01100"){
        pathview(gene.data=ranks, pathway.id=pathID, gene.idtype="kegg",
               kegg.dir=paste0(outPath,"/kegg/p",padjVal,"_lfc",lfcVal,"/",grp),
               out.suffix=paste0("GSEA_",grp),
               species="cel", limit=list(gene=max(abs(ranks)), cpd=1))
        system(paste0("mv ",outPath,"/",pathID,".GSEA_",grp,".png ",
                      outPath,"/kegg/p",padjVal,"_lfc",lfcVal,"/"))
      }
    }
  }

  mkk <- enrichMKEGG(gene = paste0("CELE_",salmonSig$sequenceID),
                     organism = 'cel',
                     qvalueCutoff = 0.1)
  #head(mkk)
  if(!is.null(mkk)){
    write.table(kk2,file=paste0(outPath,"/kegg/",fileNamePrefix,"moduleKegg_",grp,
                                "_padj", padjVal, "_lfc", lfcVal,".tsv"),sep="\t")
  }
  #browseKEGG(kk, 'cel03040')
}


# #########################-
# # wikipathway enrichment -----
# #########################-
# grp=useContrasts[1]
# for (grp in useContrasts){
#
#   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
#
#   if(filterData){
#     # remove filtered genes
#     idx<-salmon$wormbaseID %in% toFilter
#     salmon<-salmon[!idx,]
#   }
#
#   salmonSig<-getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
#                                  namePadjCol="padj",
#                                  nameLfcCol="log2FoldChange",
#                                  direction="both",
#                                  chr="all", nameChrCol="chr")
#   ranks<-salmonSig$log2FoldChange
#   names(ranks)<-salmonSig$entrezID
#   ranks<-ranks[!is.na(names(ranks))]
#   ranks<-sort(ranks,decreasing=T)
#   ewp<-enrichWP(names(ranks), organism = "Caenorhabditis elegans")
#   head(ewp)
#   dotplot(ewp)
#   upsetplot(ewp)
#   gewp<-gseWP(ranks, organism = "Caenorhabditis elegans")
#   head(gewp)
# }


#########################-
# reactomeDB enrichment -----
#########################-
library(ReactomePA)

grp=useContrasts[3]
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
  ranks<-salmonSig$log2FoldChange
  names(ranks)<-salmonSig$entrezID
  ranks<-ranks[!is.na(names(ranks))]
  ranks<-sort(ranks,decreasing=T)

  x <- enrichPathway(gene=names(ranks), organism="celegans",
                   qvalueCutoff = 0.1, readable=TRUE)
  head(x)
  plotList<-list()
  if(dim(x)[1]>0){
    plotList[["dotplot"]]<-dotplot(x, showCategory=20) + ggtitle(grp)
    plotList[["cnetstar"]] <- cnetplot(x, foldChange=ranks, colorEdge = TRUE,
                                       showCategory=10)
    plotList[["cnetcircle"]] <- cnetplot(x, foldChange=ranks, circular = TRUE,
                                         colorEdge = TRUE, showCategory=10)
    plotList[["heatplot"]] <- heatplot(x, foldChange=ranks, showCategory=30)

    x2 <- pairwise_termsim(x)
    if(dim(x2)[1]>6){
      plotList[["treeplot"]] <- treeplot(x2, hclust_method="average", fontsize=5)
    }
    plotList[["emapplot"]] <- emapplot(x2, cex_category=1,
                                       cex_label_category=0.6, cex_line=0.5) +
            ggtitle(grp)
    plotList[["upsetplot"]] <- upsetplot(x, n=10)

    pdf(file=paste0(outPath, "/kegg/",fileNamePrefix,"reactomePA_ORA_", grp,".pdf"),
        width=11,height=8,paper="a4r")
    #p<-ggpubr::ggarrange(plotlist=gseaList,nrow=2,ncol=1)
    p<-gridExtra::marrangeGrob(grobs=plotList, nrow=1, ncol=1)
    print(p)
    dev.off()

  }


  y <- gsePathway(ranks, organism="celegans",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE)
  if(dim(y)[1]>0){
    head(y)
    plotList[["dotplot"]]<-dotplot(y, showCategory=20) + ggtitle(grp)
    plotList[["cnetstar"]] <- cnetplot(y, foldChange=ranks, colorEdge = TRUE,
                                       showCategory=10)
    plotList[["cnetcircle"]] <- cnetplot(y, foldChange=ranks, circular = TRUE,
                                         colorEdge = TRUE,showCategory=10)
    plotList[["heatplot"]] <- heatplot(y, foldChange=ranks, showCategory=30)

    y2 <- pairwise_termsim(y)
    if(dim(y2)[1]>6){
      plotList[["treeplot"]] <- treeplot(y2, hclust_method = "average",fontsize=5)
    }
    plotList[["emapplot"]] <- emapplot(y2, cex_category=1,
                                       cex_label_category=0.6, cex_line=0.5) +
      ggtitle(grp)
    plotList[["upsetplot"]] <- upsetplot(y, n=10)

    pdf(file=paste0(outPath, "/kegg/",fileNamePrefix,"reactomePA_GSEA_", grp,".pdf"),
        width=11,height=8,paper="a4r")
    #p<-ggpubr::ggarrange(plotlist=gseaList,nrow=2,ncol=1)
    p<-gridExtra::marrangeGrob(grobs=plotList, nrow=1, ncol=1)
    print(p)
    dev.off()
  }
}




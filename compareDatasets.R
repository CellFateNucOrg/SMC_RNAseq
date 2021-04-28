library(readxl)
library(ggplot2)
library(EnhancedVolcano)
library(dplyr)
library(tidyr)
library(eulerr)
library(lattice)
library(gridExtra)

source("functions.R")
source("./variableSettings.R")

if(filterData){
  fileNamePrefix<-filterPrefix
}

eulerLabelsType<-c("counts","percent")
eulerLabelsType<-c("counts")

###########################
## compare samples
##########################


#######-
## venn diagrams------
#######-


## significantly changed genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))
  sigTables[[prettyGeneName(grp)]]<-as.data.frame(getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                 namePadjCol="padj",
                                 nameLfcCol="log2FoldChange",
                                 direction="both",
                                 chr="all", nameChrCol="chr"))
}

# check if datasets have chrX genes included
includeChrX<-"chrX" %in% unlist(lapply(sigTables,"[","chr"))

pdf(file=paste0(outPath, "/plots/",fileNamePrefix,"venn_allGenes_",
                paste(groupsOI, collapse="_"),"_padj",
                padjVal, "_lfc", lfcVal,".pdf"),
    width=5,height=10,paper="a4")


sigGenes<-lapply(sigTables,"[[","wormbaseID")
fit<-euler(sigGenes)
p1<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste0("All genes: |lfc|>", lfcVal, ", padj<",padjVal,"\n",
          paste(lapply(row.names(fit$ellipses), function(x){
          paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
          }), collapse="  ")), fontsize=8))

print(p1)

# dotplot(resid(fit), xlab = "Residuals",
#         panel = function(...) {
#           panel.abline(v = 0, lty = 2)
#           panel.dotplot(...)
#         })
# error_plot(fit)
# coef(fit)
#
# #http://www.pangloss.com/wiki/VennSignificance
# sigResult$phyper<-1-phyper(q=sum(subTbl[,1]*subTbl[,2]),#overlap
# m=sum(subTbl[,1]), #changed in first dataset
# n=nrow(subTbl)-sum(subTbl[,1]), #not changed in first dataset
# k=sum(subTbl[,2])) #changed in second dataset

####### chrX----
if(includeChrX){
  xchr<-lapply(sigTables,function(x) x[x$chr=="chrX",])
  sigGenes<-lapply(xchr,"[[","wormbaseID")
  fit<-euler(sigGenes)
  fit
  p2<-plot(fit, quantities=list(type=eulerLabelsType),
           main=list(label=paste0("chrX: |lfc|>", lfcVal, ", padj<",padjVal,"\n",
                       paste(lapply(row.names(fit$ellipses), function(x){
             paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
           }), collapse="  ")), fontsize=8))
} else {
  p2<-NULL
}
print(p2)

##### autosomes -----
achr<-lapply(sigTables,function(x) x[x$chr!="chrX",])
sigGenes<-lapply(achr, "[[", "wormbaseID")
fit<-euler(sigGenes)
p3<-plot(fit, quantities=list(type=eulerLabelsType),
       main=list(label=paste0("Autosomal: |lfc|>", lfcVal, ", padj<",padjVal,"\n",
            paste(lapply(row.names(fit$ellipses), function(x){
            paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
            }), collapse="  ")), fontsize=8))

print(p3)
dev.off()
#p<-arrangeGrob(grobs=list(p1,p2,p3), ncol=3,padding=2)
# p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
# ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_allGenes_",
#                                 paste(groupsOI, collapse="_"),"_padj",
#                                 padjVal, "_lfc", lfcVal,".pdf"),
#                 plot=p, device="pdf",width=29,height=16,units="cm")



## upregulated genes -----
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))
  sigTables[[prettyGeneName(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="gt",
                        chr="all", nameChrCol="chr"))
}

pdf(file=paste0(outPath, "/plots/",fileNamePrefix,"venn_upGenes_",
                paste(groupsOI,collapse="_"),"_padj",
                padjVal, "_lfc", lfcVal,".pdf"),
    width=5,height=10,paper="a4")

sigGenes<-lapply(sigTables, "[[" ,"wormbaseID")
fit<-euler(sigGenes)
p1<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste0("All genes up: lfc>", lfcVal, ", padj<",padjVal,"\n",
                     paste(lapply(row.names(fit$ellipses), function(x){
                       paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                     }), collapse="  ")), fontsize=8))
print(p1)


if(includeChrX){
  xchr<-lapply(sigTables,function(x) x[x$chr=="chrX",])
  sigGenes<-lapply(xchr, "[[","wormbaseID")
  fit<-euler(sigGenes)
  p2<-plot(fit, quantities=list(type=eulerLabelsType),
           main=list(label=paste("chrX up: lfc>", lfcVal, ", padj<",padjVal,"\n",
                       paste(lapply(row.names(fit$ellipses), function(x){
                         paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                       }), collapse="  ")), fontsize=8))

} else {
  p2<-NULL
}
print(p2)

achr<-lapply(sigTables,function(x) x[x$chr!="chrX",])
sigGenes<-lapply(achr, "[[","wormbaseID")
fit<-euler(sigGenes)
p3<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste("Autosomal up: lfc>", lfcVal, ", padj<",padjVal,"\n",
                     paste(lapply(row.names(fit$ellipses), function(x){
                       paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                     }), collapse="  ")), fontsize=8))
print(p3)
dev.off()

# p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
# ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_upGenes_",
#                                 paste(groupsOI,collapse="_"),"_padj",
#                                 padjVal, "_lfc", lfcVal,".pdf"),
#                 plot=p, device="pdf",width=29,height=16,units="cm")



## downregulated genes -----
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="lt",
                        chr="all", nameChrCol="chr"))
}


pdf(file=paste0(outPath, "/plots/",fileNamePrefix,"venn_downGenes_",
                paste(groupsOI,collapse="_"),"_padj",
                padjVal, "_lfc", lfcVal,".pdf"),
    width=5,height=10,paper="a4")

sigGenes<-lapply(sigTables, "[[","wormbaseID")
fit<-euler(sigGenes)
p1<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste("All genes down: lfc< -", lfcVal, ", padj<",padjVal,"\n",
                     paste(lapply(row.names(fit$ellipses), function(x){
                       paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                     }), collapse="  ")), fontsize=8))
print(p1)


if(includeChrX){
  xchr<-lapply(sigTables,function(x) x[x$chr=="chrX",])
  sigGenes<-lapply(xchr, "[[","wormbaseID")
  fit<-euler(sigGenes)
  p2<-plot(fit, quantities=list(type=eulerLabelsType),
           main=list(label=paste("chrX down: lfc< -", lfcVal, ", padj<",padjVal,"\n",
                       paste(lapply(row.names(fit$ellipses), function(x){
                         paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                       }), collapse="  ")), fontsize=8))

} else {
  p2<-NULL
}
print(p2)

achr<-lapply(sigTables,function(x) x[x$chr!="chrX",])
sigGenes<-lapply(achr, "[[","wormbaseID")
fit<-euler(sigGenes)
p3<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste("Autosomal down: lfc< -", lfcVal, ", padj<",padjVal,"\n",
                     paste(lapply(row.names(fit$ellipses), function(x){
                       paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                     }), collapse="  ")), fontsize=8))
print(p3)
dev.off()

# p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
# ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_downGenes_",
#                                 paste(groupsOI,collapse="_"),"_padj",
#                                 padjVal, "_lfc", lfcVal,".pdf"),
#                 plot=p, device="pdf",width=29,height=16,units="cm")




########-
## significantly changed count - barplot ------
########-
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(getSignificantGenes(salmon,
                                                  padj=padjVal, lfc=lfcVal,
                                                  namePadjCol="padj",
                                                  nameLfcCol="log2FoldChange",
                                                  direction="both",
                                                  chr="all", nameChrCol="chr"))
}

sigPerChr<-lapply(sigTables, "[", ,"chr")
sigPerChr<-as.data.frame(do.call(rbind,lapply(sigPerChr,table)))
sigPerChr$SMC<-rownames(sigPerChr)
sigPerChr<-sigPerChr %>% gather(colnames(sigPerChr)[1:6],key=chr, value=genes)
sigPerChr$chr<-gsub("chr","",sigPerChr$chr)


p1<-ggplot(sigPerChr, aes(x=chr,y=genes,group=SMC)) + facet_grid(cols=vars(SMC)) +
  geom_bar(stat="identity",position=position_dodge(),aes(fill=chr), show.legend = FALSE) +
  geom_text(aes(label=genes), vjust=-0.2, color="black",
          position = position_dodge(0.9), size=3) + ylab("Number of genes") +
  ggtitle("Number of significantly changed genes per chromosome")  +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2)


########-
## up regulated count - barplot-----
########-

sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(getSignificantGenes(salmon,
                                                      padj=padjVal, lfc=lfcVal,
                                                      namePadjCol="padj",
                                                      nameLfcCol="log2FoldChange",
                                                      direction="gt",
                                                      chr="all",
                                                      nameChrCol="chr"))
}


sigPerChr<-lapply(sigTables, "[", ,"chr")
sigPerChr<-as.data.frame(do.call(rbind,lapply(sigPerChr,table)))
sigPerChr$SMC<-rownames(sigPerChr)
sigPerChr<-sigPerChr %>% gather(colnames(sigPerChr)[1:6],key=chr, value=genes)
sigPerChr$chr<-gsub("chr","",sigPerChr$chr)
ymax<-max(sigPerChr$genes)

p2<-ggplot(sigPerChr, aes(x=chr,y=genes,group=SMC)) + facet_grid(cols=vars(SMC)) +
  geom_bar(stat="identity",position=position_dodge(),aes(fill=chr), show.legend = FALSE) +
  geom_text(aes(label=genes), vjust=-0.2, color="black",
            position = position_dodge(0.9), size=3) + ylab("Number of genes") +
  ggtitle("Number of significantly upregulated genes per chromosome")  +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2)

########-
## down regulated count -barplot-----
########-

sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(getSignificantGenes(
                                        salmon, padj=padjVal,
                                        lfc= lfcVal,
                                        namePadjCol="padj",
                                        nameLfcCol="log2FoldChange",
                                        direction="lt",
                                        chr="all",
                                        nameChrCol="chr"))
}

sigPerChr<-lapply(sigTables, "[", ,"chr")
sigPerChr<-as.data.frame(do.call(rbind,lapply(sigPerChr,table)))
sigPerChr$SMC<-rownames(sigPerChr)
sigPerChr<-sigPerChr %>% gather(colnames(sigPerChr)[1:6],key=chr, value=genes)
sigPerChr$chr<-gsub("chr","",sigPerChr$chr)

p3<-ggplot(sigPerChr, aes(x=chr,y=genes,group=SMC)) +
  facet_grid(cols=vars(SMC), switch="x") +
  geom_bar(stat="identity",position=position_dodge(),aes(fill=chr), show.legend = FALSE) +
  geom_text(aes(label=genes), vjust=1.2, color="black",
            position = position_dodge(0.9), size=3) + ylab("Number of genes") +
  ggtitle("Number of significantly downregulated genes per chromosome")  +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2) +
  scale_y_reverse(limits=c(ymax,0)) + scale_x_discrete(position = "top")

p<-ggpubr::ggarrange(p1,p2,p3,ncol=1,nrow=3)
ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,
                                "bar_countsPerChr_",paste(groupsOI,
                                collapse="_"),"_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=19,height=29,units="cm")



########-
## significantly changed LFC &count- boxplot&barplot-----
########-
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds"))
  salmon<-salmon[!is.na(salmon$chr),]
  rownames(salmon)<-NULL
  sigTables[[grp]]<-as.data.frame(getSignificantGenes(salmon,
                                        padj=padjVal, lfc=lfcVal,
                                        namePadjCol="padj",
                                        nameLfcCol="log2FoldChange",
                                        direction="both",
                                        chr="all", nameChrCol="chr"))
}



# sigPerChr<-lapply(sigTables, "[", ,c("wormbaseID","chr","log2FoldChange"))
# for(g in names(sigPerChr)){ sigPerChr[[g]]$SMC<-g }
# sigPerChr<-as.data.frame(do.call(rbind,sigPerChr))
# #sigPerChr$SMC<-gsub("\\.\\d*$","",rownames(sigPerChr))
# #rownames(sigPerChr)<-NULL
# sigPerChr$updown<-"both"
#
# yminmax=c(0,median(abs(sigTbl$log2FoldChange))+quantile(abs(sigTbl$log2FoldChange))[4]*2)
# p1<-ggplot(sigPerChr,aes(x=chr,y=log2FoldChange,fill=SMC))+
#   geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2),outlier.shape=NA)+ ylim(yminmax) +
#   ggtitle("LFC of significantly changed genes by chromosome") +
#   theme_minimal() + scale_fill_brewer(palette="Dark2")


# upregulated
sigList<-lapply(sigTables, getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="gt")

sigList<-lapply(sigList, "[", ,c("wormbaseID","chr","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList$updown<-"up"
sigTbl<-sigList


# downregulated
sigList<-lapply(sigTables, getSignificantGenes,
                padj=padjVal,lfc=lfcVal,direction="lt")

sigList<-lapply(sigList, "[", ,c("wormbaseID","chr","log2FoldChange"))
for(g in names(sigList)){ sigList[[g]]$SMC<-g }
sigList<-do.call(rbind,sigList)
sigList$updown<-"down"
sigTbl<-rbind(sigTbl,sigList)
sigTbl$chr<-as.factor(sigTbl$chr)
sigTbl$SMC<-as.factor(sigTbl$SMC)
sigTbl$updown<-factor(sigTbl$updown, levels=c("up","down"))
rownames(sigTbl)<-NULL

yminmax=c(0,median(abs(sigTbl$log2FoldChange))+quantile(abs(sigTbl$log2FoldChange))[4]*2)
p2<-ggplot(sigTbl,aes(x=chr,y=abs(log2FoldChange),fill=SMC)) +
  geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2),outlier.shape=NA,lwd=0.1,fatten=3) +
  facet_grid(cols=vars(updown)) + ylim(yminmax) +
  ggtitle("Absolute LFC of significantly changed genes by chromosome") +
  theme_minimal() + scale_fill_brewer(palette="Dark2")

yminmax=c(0,median(abs(sigTbl$log2FoldChange))+quantile(abs(sigTbl$log2FoldChange))[4]*2)
p3<-ggplot(sigTbl,aes(x=updown,y=abs(log2FoldChange),fill=SMC)) +
  geom_boxplot(notch=T, varwidth=F, position=position_dodge2(padding=0.2), outlier.shape=NA) +
  ggtitle("LFC of significantly changed autosomal genes") +
  theme_minimal()  +
  scale_y_continuous(limits = yminmax) + xlab(NULL) +
  scale_fill_brewer(palette="Dark2")

countbychr<-sigTbl %>% group_by(updown,chr,SMC) %>%
  dplyr::summarise(count=n(),groups="keep")
yminmax=c(0,max(countbychr$count))
p4<-ggplot(countbychr,aes(x=chr,y=count,fill=SMC)) +
  geom_bar(stat="identity",position=position_dodge()) +
  facet_grid(cols=vars(updown))+ ylim(yminmax) +
  ggtitle("Count of significantly changed genes by chromosome") +
  theme_minimal()  + scale_fill_brewer(palette="Dark2")

yminmax=c(0,countbychr$count[order(countbychr$count,decreasing=T)[2]])
p5<-ggplot(countbychr,aes(x=chr,y=count,fill=SMC)) +
  geom_bar(stat="identity",position=position_dodge()) +
  facet_grid(cols=vars(updown))+
  ggtitle("Count of significantly changed genes by chromosome") +
  theme_minimal()  + scale_y_continuous(limits = yminmax) +
  scale_fill_brewer(palette="Dark2")


countbytype<-sigTbl %>% filter(chr!="chrX") %>% group_by(updown,SMC) %>%
  dplyr::summarise(count=n(),groups="keep")

yminmax=c(0,countbytype$count[order(countbytype$count,decreasing=T)[1]])
p6<-ggplot(countbytype,aes(x=updown,y=count,fill=SMC)) +
  geom_bar(stat="identity",position=position_dodge(),lwd=0.1) +
  ylim(yminmax) +
  ggtitle("Count of significantly changed autosomal genes") +
  theme_minimal() +xlab(NULL) +
  scale_fill_brewer(palette="Dark2") + ylab("Number of genes")

# up vs down
bb<-countbytype %>% pivot_wider(names_from=updown,values_from=count) %>% group_by(SMC) %>% mutate(ratioupdown=as.character(round(up/down,1)), updown="up")
bb1<-bb
bb1$updown<-"down"
bb1$ratioupdown<-"1"

# diff datasets
aa<-countbytype %>% group_by(updown) %>% mutate(ratio=round(count/min(count),1))

p6<-p6+geom_text(data=aa,aes(x=updown,y=10,label=ratio), color="black",
             position = position_dodge(1)) +
  geom_text(data=rbind(bb,bb1),aes(x=updown,y=max(bb$up),label=ratioupdown), color="black",
                    position = position_dodge(1),vjust=-0.5)
p<-ggpubr::ggarrange(p2,p4,p5,ncol=1,nrow=3)
ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,
                                "count&LFCbyChr_",
                                paste(groupsOI,collapse="_"), "_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=19,height=29,units="cm")


p<-ggpubr::ggarrange(p3,p6,ncol=1,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,
                                "count&LFCautosomes_",
                                paste(groupsOI,collapse="_"), "_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=13,height=22,units="cm")



#########-
## correlation-----
#########-

geneTable<-NULL
for (grp in groupsOI){
  salmon<-as.data.frame(readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds")))
  colnames(salmon)[colnames(salmon)=="log2FoldChange"]<-paste0(grp,"_lfc")
  if(is.null(geneTable)){
    geneTable<-as.data.frame(salmon)[,c("wormbaseID","chr",paste0(grp,"_lfc"))]
  } else {
    geneTable<-full_join(geneTable,salmon[,c("wormbaseID","chr",paste0(grp,"_lfc"))], by=c("wormbaseID","chr"))
  }
}

combnTable<-combn(1:length(groupsOI),m=2)
# all genes
# geneTable<-na.omit(geneTable)
# gt<-geneTable[geneTable$kle2cs_lfc<1 & geneTable$scc1cs_lfc>2,]
# dim(gt)
# david<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/david_wbid2entrez_WS278.txt")
# showPeter<-david[david$From %in% gt$wormbaseID,]
# write.table(showPeter,file=paste0(outPath,"/txt/klelt1_sccgt2.tsv"),sep="\t")

lfcCols<-grep("_lfc$",names(geneTable))
minScale<-min(geneTable[,lfcCols])*1.1
maxScale<-max(geneTable[,lfcCols])*1.1

if(plotPDFs==T){
  pdf(file=paste0(outPath, "/plots/",fileNamePrefix,"cor_allGenes.pdf"),
      width=5, height=5, paper="a4")
}
for (i in 1:ncol(combnTable)){
  grp1<-groupsOI[combnTable[1,i]]
  grp2<-groupsOI[combnTable[2,i]]

  if(plotPDFs==F){
    png(file=paste0(outPath, "/plots/",fileNamePrefix,"cor_allGenes_",grp1,"_",
                    grp2,".png"), width=5, height=5, units="in", res=150)
  }
  Rval<-round(cor(geneTable[,paste0(grp1,"_lfc")],
                  geneTable[,paste0(grp2,"_lfc")]),2)
  #smoothScatter(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")],
  #              xlab=grp1,ylab=grp2,xlim=c(minScale,maxScale), nrpoints=1000,
  #             col="red", colramp = colorRampPalette(c("white", rev(grey.colors(10)))),
  #              ylim=c(minScale,maxScale),transformation=function(x) x^.25)
  plot(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")],pch=16,
       cex=0.5,col="#11111155",xlab=prettyGeneName(grp1),
       ylab=prettyGeneName(grp2), xlim=c(minScale,maxScale),
       ylim=c(minScale,maxScale))
  abline(v=0,h=0,col="grey60",lty=3)
  bestFitLine<-lm(geneTable[,paste0(grp2,"_lfc")]~geneTable[,paste0(grp1,"_lfc")])
  abline(bestFitLine,col="red")
  title(main=paste0("All genes ",prettyGeneName(grp1)," vs ",
                    prettyGeneName(grp2)," (R=",Rval,")"),
        sub=paste0(nrow(geneTable)," genes"))
  if(plotPDFs==F){
    dev.off()
  }
}

tmp<-geneTable # create backup copy

###########-
## correlation of chrX genes-----
###########-
if(includeChrX){
  geneTable<-na.omit(geneTable)
  geneTable<-geneTable[geneTable$chr=="chrX",]
  for (i in 1:ncol(combnTable)){
    grp1<-groupsOI[combnTable[1,i]]
    grp2<-groupsOI[combnTable[2,i]]
    if(plotPDFs==F){
      png(file=paste0(outPath, "/plots/",fileNamePrefix,"cor_chrX_",grp1,"_",
                      grp2,".png"), width=5,
          height=5, units="in", res=150)
    }
    Rval<-round(cor(geneTable[,paste0(grp1,"_lfc")],
                    geneTable[,paste0(grp2,"_lfc")]),2)
    #smoothScatter(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")],
    #              xlab=grp1,ylab=grp2,xlim=c(minScale,maxScale),
    #              ylim=c(minScale,maxScale),transformation=function(x) x^.25)
    plot(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")],pch=16,
         cex=0.5,col="#11111155",xlab=prettyGeneName(grp1),
         ylab=prettyGeneName(grp2), xlim=c(minScale,maxScale),
         ylim=c(minScale,maxScale))
    abline(v=0,h=0,col="grey60",lty=3)
    bestFitLine<-lm(geneTable[,paste0(grp2,"_lfc")]~geneTable[,paste0(grp1,"_lfc")])
    abline(bestFitLine,col="red")
    title(main=paste0("chrX genes ",prettyGeneName(grp1)," vs ",
                 prettyGeneName(grp2)," (R=",Rval,")"),
          sub=paste0(nrow(geneTable)," genes"))
    if(plotPDFs==F){
      dev.off()
    }
  }
}


###########-
## correlation of autosomal genes-----
###########-
geneTable<-tmp # recover original geneTable
geneTable<-geneTable[geneTable$chr!="chrX",]
for (i in 1:ncol(combnTable)){
  grp1<-groupsOI[combnTable[1,i]]
  grp2<-groupsOI[combnTable[2,i]]
  if(plotPDFs==F){
    png(file=paste0(outPath, "/plots/",fileNamePrefix,"cor_autosomal_",grp1,"_",
                    grp2,".png"), width=5,height=5, units="in", res=150)
  }
  Rval<-round(cor(geneTable[,paste0(grp1,"_lfc")],
                  geneTable[,paste0(grp2,"_lfc")]),2)
  plot(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")],pch=16,
       cex=0.5, col="#11111155",xlab=prettyGeneName(grp1),
       ylab=prettyGeneName(grp2), xlim=c(minScale,maxScale),
       ylim=c(minScale,maxScale))
  #smoothScatter(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")],
  #              xlab=grp1,ylab=grp2,xlim=c(minScale,maxScale),
  #              ylim=c(minScale,maxScale),transformation=function(x) x^.25)
  abline(v=0,h=0,col="grey60",lty=3)
  bestFitLine<-lm(geneTable[,paste0(grp2,"_lfc")]~geneTable[,paste0(grp1,"_lfc")])
  abline(bestFitLine,col="red")
  title(main=paste0("Autosomal genes ",prettyGeneName(grp1)," vs ",
               prettyGeneName(grp2)," (R=",Rval,")"),
        sub=paste0(nrow(geneTable)," genes"))
  if(plotPDFs==F){
    dev.off()
  }
}

if(plotPDFs==T){
  dev.off()
}



if(!combineChrAX){
  #########################
  ## compare LFC to wt mean
  #########################
  dds<-readRDS(file=paste0(outPath,"/rds/dds_object.rds"))
  wtMean<-rowMeans(counts(dds)[,colData(dds)$SMC==controlGrp])
  idx<-wtMean!=0
  logwtcounts<-log2(wtMean[idx])
  logwtcounts[is.infinite(logwtcounts)]<-NA

  geneTable<-NULL
  for (grp in groupsOI){
    salmon<-as.data.frame(readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults_p",padjVal,".rds")))
    colnames(salmon)[colnames(salmon)=="log2FoldChange"]<-paste0(grp,"_lfc")
    if(is.null(geneTable)){
      geneTable<-as.data.frame(salmon)[,c("wormbaseID","chr", paste0(grp,"_lfc"))]
    } else {
      geneTable<-inner_join(geneTable,salmon[,c("wormbaseID", "chr",paste0(grp,"_lfc"))], by=c("wormbaseID","chr"))
    }
  }

  dim(geneTable)
  # all genes
  geneTable<-na.omit(geneTable)

  lfcCols<-grep("_lfc$",names(geneTable))
  minScale<-min(geneTable[,lfcCols])*1.1
  maxScale<-max(geneTable[,lfcCols])*1.1

  if(plotPDFs==T){
    pdf(file=paste0(outPath, "/plots/",fileNamePrefix,"cor_wtExpr.pdf"), width=5,
        height=5, paper="a4")
  }
  for (grp in groupsOI){
    if(plotPDFs==F){
      png(file=paste0(outPath, "/plots/",fileNamePrefix,"cor_wtExpr_",grp1,"_",
                      grp2,".png"), width=5, height=5, units="in", res=150)
    }
    idx<-match(geneTable$wormbaseID,names(wtMean))
    noNAtbl<-na.omit(cbind(as.vector(geneTable[,paste0(grp,"_lfc")]),
                           logwtcounts[idx]))
    Rval<-round(cor(noNAtbl[,1],noNAtbl[,2]),2)

    plot(logwtcounts[idx],geneTable[,paste0(grp,"_lfc")],pch=16,
         cex=0.5,col="#11111155", xlim=c(minScale,maxScale),
         ylab=paste0("log2(",prettyGeneName(grp),"/", controlGrp,")"),
         xlab=paste0("log2(",controlGrp,")"))
    abline(v=0,h=0,col="grey60",lty=3)
    bestFitLine<-lm(geneTable[,paste0(grp,"_lfc")]~logwtcounts[idx],
                    na.action=na.exclude)
    abline(bestFitLine,col="red")
    title(paste0(prettyGeneName(grp)," LFC vs ",controlGrp,
                 " log counts (R=", Rval,")"))
    if(plotPDFs==F){
      dev.off()
    }
  }
  dev.off()
}


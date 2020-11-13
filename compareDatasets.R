library(readxl)
library(ggVennDiagram)
library(ggplot2)
library(EnhancedVolcano)

source("functions.R")

outPath="."
padjVal=0.05
lfcVal=0.5

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
  salmon<-readRDS(paste0(outPath,"/rds/salmon_",grp,"_DESeq2_fullResults.rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                 namePadjCol="padj",
                                 nameLfcCol="log2FoldChange",
                                 direction="both",
                                 chr="all", nameChrCol="chr"))
}

sigGenes<-lapply(sigTables, "[", ,"wormbaseID")
p1<-ggVennDiagram(sigGenes) + ggtitle(label=paste0("All genes: |lfc|>", lfcVal, ", padj<",padjVal))

xchr<-lapply(sigTables,function(x) x[x$chr=="chrX",])
sigGenes<-lapply(xchr, "[", ,"wormbaseID")
p2<-ggVennDiagram(sigGenes) + ggtitle(label=paste0("chrX genes: |lfc|>", lfcVal, ", padj<",padjVal))

achr<-lapply(sigTables,function(x) x[x$chr!="chrX",])
sigGenes<-lapply(achr, "[", ,"wormbaseID")
p3<-ggVennDiagram(sigGenes) + ggtitle(label=paste0("Autosomal genes: |lfc|>", lfcVal, ", padj<",padjVal))

p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_",paste(groupsOI,
                                                  collapse="_"),"_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-1",".pdf"),
                plot=p, device="pdf",width=29,height=16,units="cm")


########
## significantly changed count - barplot
########
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


########
## up regulated count - barplot
########

sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/salmon_",grp,"_DESeq2_fullResults.rds"))

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

p2<-ggplot(sigPerChr, aes(x=chr,y=genes,group=SMC)) + facet_grid(cols=vars(SMC)) +
  geom_bar(stat="identity",position=position_dodge(),aes(fill=chr), show.legend = FALSE) +
  geom_text(aes(label=genes), vjust=-0.2, color="black",
            position = position_dodge(0.9), size=3) + ylab("Number of genes") +
  ggtitle("Number of significantly upregulated genes per chromosome")  +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2)

########
## down regulated count -barplot
########

sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/salmon_",grp,"_DESeq2_fullResults.rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(getSignificantGenes(salmon,
                                                                      padj=padjVal, lfc= -lfcVal,
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

p3<-ggplot(sigPerChr, aes(x=chr,y=genes,group=SMC)) + facet_grid(cols=vars(SMC)) +
  geom_bar(stat="identity",position=position_dodge(),aes(fill=chr), show.legend = FALSE) +
  geom_text(aes(label=genes), vjust=-0.2, color="black",
            position = position_dodge(0.9), size=3) + ylab("Number of genes") +
  ggtitle("Number of significantly downregulated genes per chromosome")  +
  theme_minimal() + scale_fill_grey(start=0.8, end=0.2)

p<-ggpubr::ggarrange(p1,p2,p3,ncol=1,nrow=3)
ggplot2::ggsave(filename=paste0(outPath, "/plots/bar_countsPerChr_",paste(groupsOI,
                                collapse="_"),"_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=19,height=29,units="cm")



#########
## correlation
#########

geneTable<-NULL
for (grp in groupsOI){
  salmon<-as.data.frame(readRDS(paste0(outPath,"/rds/salmon_",grp,"_DESeq2_fullResults.rds")))
  colnames(salmon)[colnames(salmon)=="log2FoldChange"]<-paste0(grp,"_lfc")
  if(is.null(geneTable)){
    geneTable<-as.data.frame(salmon)[,c("wormbaseID","chr",paste0(grp,"_lfc"))]
  } else {
    geneTable<-inner_join(geneTable,salmon[,c("wormbaseID","chr",paste0(grp,"_lfc"))], by=c("wormbaseID","chr"))
  }
}

combnTable<-combn(1:length(groupsOI),m=2)
# all genes
geneTable<-na.omit(geneTable)

for (i in 1:ncol(combnTable)){
  grp1<-groupsOI[combnTable[1,i]]
  grp2<-groupsOI[combnTable[2,i]]
  png(file=paste0(outPath, "/plots/cor_allGenes_",grp1,"_",grp2,".png"), width=5,
      height=5, units="in", res=150)
  minScale<-min(geneTable[,c(paste0(grp1,"_lfc"),paste0(grp2,"_lfc"))])
  maxScale<-max(geneTable[,c(paste0(grp1,"_lfc"),paste0(grp2,"_lfc"))])
  Rval<-round(cor(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")]),2)
  plot(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")],pch=1,
       cex=0.5,col="#111111ee",xlab=grp1,ylab=grp2,xlim=c(minScale,maxScale),
       ylim=c(minScale,maxScale))
  bestFitLine<-lm(geneTable[,paste0(grp2,"_lfc")]~geneTable[,paste0(grp1,"_lfc")])
  abline(bestFitLine,col="red")
  title(paste0("All genes ",grp1," vs ", grp2," (R=",Rval,")"))
  dev.off()
}

geneTable<-na.omit(geneTable)
tmp<-geneTable
geneTable<-geneTable[geneTable$chr=="chrX",]
for (i in 1:ncol(combnTable)){
  grp1<-groupsOI[combnTable[1,i]]
  grp2<-groupsOI[combnTable[2,i]]
  png(file=paste0(outPath, "/plots/cor_chrX_",grp1,"_",grp2,".png"), width=5,
      height=5, units="in", res=150)
  minScale<-min(geneTable[,c(paste0(grp1,"_lfc"),paste0(grp2,"_lfc"))])
  maxScale<-max(geneTable[,c(paste0(grp1,"_lfc"),paste0(grp2,"_lfc"))])
  Rval<-round(cor(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")]),2)
  plot(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")],pch=1,
       cex=0.5,col="#111111ee",xlab=grp1,ylab=grp2,xlim=c(minScale,maxScale),
       ylim=c(minScale,maxScale))
  bestFitLine<-lm(geneTable[,paste0(grp2,"_lfc")]~geneTable[,paste0(grp1,"_lfc")])
  abline(bestFitLine,col="red")
  title(paste0("chrX genes ",grp1," vs ", grp2," (R=",Rval,")"))
  dev.off()
}

geneTable<-tmp
geneTable<-geneTable[geneTable$chr!="chrX",]
for (i in 1:ncol(combnTable)){
  grp1<-groupsOI[combnTable[1,i]]
  grp2<-groupsOI[combnTable[2,i]]
  png(file=paste0(outPath, "/plots/cor_autosomal_",grp1,"_",grp2,".png"), width=5,
      height=5, units="in", res=150)
  minScale<-min(geneTable[,c(paste0(grp1,"_lfc"),paste0(grp2,"_lfc"))])
  maxScale<-max(geneTable[,c(paste0(grp1,"_lfc"),paste0(grp2,"_lfc"))])
  Rval<-round(cor(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")]),2)
  plot(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")],pch=1,
       cex=0.5,col="#111111ee",xlab=grp1,ylab=grp2,xlim=c(minScale,maxScale),
       ylim=c(minScale,maxScale))
  bestFitLine<-lm(geneTable[,paste0(grp2,"_lfc")]~geneTable[,paste0(grp1,"_lfc")])
  abline(bestFitLine,col="red")
  title(paste0("Autosomal genes ",grp1," vs ", grp2," (R=",Rval,")"))
  dev.off()
}




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

if(filterData){
  fileNamePrefix=filterPrefix
}

#####################################################
## compare to germline -soma data
#####################################################

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


for (grp in groupsOI){
  #grp="dpy26cs"
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

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


  x<-list(salmon=salmonSig$wormbaseID, germline=reinke$wormbaseID)
  names(x)<-c(prettyGeneName(grp),"germline")
  p1<-ggVennDiagram(x) + ggtitle(label=paste0(grp," vs Reinke(2004): |lfc|>", lfcVal, ", padj<",padjVal))

  x<-list(salmon=salmonSig$wormbaseID,
          germlineL4=boeck$wormbaseID[boeck$germline=="germlineL4"],
          somaL4=boeck$wormbaseID[boeck$germline=="somaL4"])
  names(x)<-c(prettyGeneName(grp),"germlineL4","somaL4")
  txtLabels<-list()
  txtLabels[paste0("% ",names(x)[1]," in ",names(x)[2])]<-round(100*length(intersect(x[[1]],x[[2]]))/length(x[[1]]),1)
  txtLabels[paste0("% ",names(x)[1]," in ",names(x)[3])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[1]]),1)
  p2<-ggVennDiagram(x) + ggtitle(label=paste0(grp," vs Boeck(2016): |lfc|>", lfcVal, ", padj<",padjVal), subtitle=paste0(txtLabels[[1]], names(txtLabels)[1], " & \n", txtLabels[[2]], names(txtLabels)[2]))

  x<-list(salmon=salmonSig$wormbaseID,
          gonadYA=boeck$wormbaseID[boeck$germline=="gonadYA"],
          somaYA=boeck$wormbaseID[boeck$germline=="somaYA"])
  names(x)<-c(prettyGeneName(grp),"gonadYA","somaYA")
  txtLabels<-list()
  txtLabels[paste0("% ",names(x)[1]," in ",names(x)[2])]<-round(100*length(intersect(x[[1]],x[[2]]))/length(x[[1]]),1)
  txtLabels[paste0("% ",names(x)[1]," in ",names(x)[3])]<-round(100*length(intersect(x[[1]],x[[3]]))/length(x[[1]]),1)
  p3<-ggVennDiagram(x) + ggtitle(label=paste0(grp," vs Boeck(2016): |lfc|>", lfcVal, ", padj<",padjVal), subtitle=paste0(txtLabels[[1]], names(txtLabels)[1], " & \n", txtLabels[[2]], names(txtLabels)[2]))

  p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_",
                                  grp,"VsGermline_padj",
                                  padjVal,"_lfc", lfcVal,".pdf"),
                  plot=p, device="pdf",width=29,height=11,units="cm")
}



## upregulated genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="gt",
                        chr="all", nameChrCol="chr"))
}

sigGenes<-lapply(sigTables, "[", ,"wormbaseID")
sigGenes[["germline"]]<-reinke$wormbaseID
p1<-ggVennDiagram(sigGenes[c(1,4)]) + ggtitle(label=paste0(grp, " genes up: lfc>", lfcVal, ", padj<",padjVal))
p2<-ggVennDiagram(sigGenes[c(2,4)]) + ggtitle(label=paste0(grp," genes up: lfc>", lfcVal, ", padj<",padjVal))
p3<-ggVennDiagram(sigGenes[c(3,4)]) + ggtitle(label=paste0(grp," genes up: lfc>", lfcVal, ", padj<",padjVal))

p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,
                                "venn_UpVsGermline_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=29,height=11,units="cm")



## upregulated genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

  sigTables[[prettyGeneName(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="lt",
                        chr="all", nameChrCol="chr"))
}

sigGenes<-lapply(sigTables, "[", ,"wormbaseID")
sigGenes[["germline"]]<-reinke$wormbaseID
p1<-ggVennDiagram(sigGenes[c(1,4)]) + ggtitle(label=paste0(grp, " genes down: lfc < -", lfcVal, ", padj<",padjVal))
p2<-ggVennDiagram(sigGenes[c(2,4)]) + ggtitle(label=paste0(grp," genes down: lfc < -", lfcVal, ", padj<",padjVal))
p3<-ggVennDiagram(sigGenes[c(3,4)]) + ggtitle(label=paste0(grp," genes down: lfc < -", lfcVal, ", padj<",padjVal))

p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)

ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,
                                "venn_DownVsGermline_padj",
                                padjVal,"_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=29,height=11,units="cm")

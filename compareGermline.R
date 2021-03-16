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


#######
## venn diagrams
#######

sigTables<-list()
sigGR<-list()
df<-data.frame(matrix(nrow=10,ncol=3))
names(df)<-groupsOI
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

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
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

  if(filterData){
    # remove filtered genes
    idx<-salmon$wormbaseID %in% toFilter
    salmon<-salmon[!idx,]

    fileNamePrefix=filterPrefix
  }

  #####################################################
  ## compare to germline -soma data
  #####################################################

  boeck<-read.csv(paste0(outPath,"/publicData/germlineSomaGenes_Boeck2016.csv"),
                   stringsAsFactors=F)
  reinke<-read.csv(paste0(outPath,"/publicData/germlineGenes_Reinke2004.csv"),
                   stringsAsFactors=F)

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
  p2<-ggVennDiagram(x) + ggtitle(label=paste0(grp," vs Boeck(2016): |lfc|>", lfcVal, ", padj<",padjVal))

  x<-list(salmon=salmonSig$wormbaseID,
          gonadYA=boeck$wormbaseID[boeck$germline=="gonadYA"],
          somaYA=boeck$wormbaseID[boeck$germline=="somaYA"])
  names(x)<-c(prettyGeneName(grp),"gonadYA","somaYA")
  p3<-ggVennDiagram(x) + ggtitle(label=paste0(grp," vs Boeck(2016): |lfc|>", lfcVal, ", padj<",padjVal))

  p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_",
                                  grp,"VsGermline_padj",
                                  padjVal,"_lfc", lfcVal,".pdf"),
                  plot=p, device="pdf",width=29,height=11,units="cm")
}

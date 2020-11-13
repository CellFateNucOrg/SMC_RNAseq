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
p1<-ggVennDiagram(sigGenes) + ggtitle(label=paste0(paste(names(sigGenes),collapse=" vs "),": |lfc|>", lfcVal, ", padj<",padjVal))

xchr<-lapply(sigTables,function(x) x[x$chr=="chrX",])
sigGenes<-lapply(xchr, "[", ,"wormbaseID")
p2<-ggVennDiagram(sigGenes) + ggtitle(label=paste0("chrX genes: |lfc|>", lfcVal, ", padj<",padjVal))

achr<-lapply(sigTables,function(x) x[x$chr!="chrX",])
sigGenes<-lapply(achr, "[", ,"wormbaseID")
p3<-ggVennDiagram(sigGenes) + ggtitle(label=paste0("autosomal genes: |lfc|>", lfcVal, ", padj<",padjVal))

p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_",paste(groupsOI,collapse="vs"),"_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-1",".pdf"),
                plot=p, device="pdf",width=19,height=20,units="cm")


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
ggplot2::ggsave(filename=paste0(outPath, "/plots/bar_PerChr",paste(groupsOI,
                                collapse="vs"),"_padj",
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
geneTable<-na.omit(geneTable)

for (i in 1:ncol(combnTable)){
  grp1<-groupsOI[combnTable[1,i]]
  grp2<-groupsOI[combnTable[2,i]]
  png(file=paste0(outPath, "/plots/cor_",grp1,"_",grp2,".png"), width=5,
      height=5, units="in", res=150)
  minScale<-min(geneTable[,c(paste0(grp1,"_lfc"),paste0(grp2,"_lfc"))])
  maxScale<-max(geneTable[,c(paste0(grp1,"_lfc"),paste0(grp2,"_lfc"))])
  Rval<-round(cor(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")]),2)
  plot(geneTable[,paste0(grp1,"_lfc")],geneTable[,paste0(grp2,"_lfc")],pch=1,
       cex=0.5,col="#111111ee",xlab=grp1,ylab=grp2,xlim=c(minScale,maxScale),
       ylim=c(minScale,maxScale))
  bestFitLine<-lm(geneTable[,paste0(grp2,"_lfc")]~geneTable[,paste0(grp1,"_lfc")])
  abline(bestFitLine,col="red")
  title(paste0(grp1," vs ", grp2," (R=",Rval,")"))
  dev.off()
}




for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/salmon_",grp,"_DESeq2_fullResults.rds"))

  #####################################################
  ## compare to public data
  #####################################################

  kramerURL<-"https://doi.org/10.1371/journal.pgen.1005698.s011"
  kramerFileName="Kramer_2015_PlotGen_S3file.xlsx"
  if(!file.exists(kramerFileName)){
    download.file(url=kramerURL,destfile=kramerFileName)
  }
  kramer<-read_excel(kramerFileName,col_types=c(rep("text",3),rep("numeric",30)))
  names(kramer)
  kramer<-kramer[,c(1:3,grep("_L3_",names(kramer)))]

  lfcVal<-0
  padjVal<-0.05

  kramerDpy27<-getSignificantGenes(kramer, padj=padjVal, lfc=lfcVal,
                      namePadjCol="dpy27_RNAi_L3_padj",
                      nameLfcCol="dpy27_RNAi_L3_log2_fold_change",
                      direction="both",
                      chr="all", nameChrCol="chr", outPath=".")

  kramerDpy21<-getSignificantGenes(kramer, padj=padjVal, lfc=lfcVal,
                                   namePadjCol="dpy21_mutant_L3_padj",
                                   nameLfcCol="dpy21_mutant_L3_log2_fold_change",
                                   direction="both",
                                   chr="all", nameChrCol="chr")

  salmonSig<-getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                 namePadjCol="padj",
                                 nameLfcCol="log2FoldChange",
                                 direction="both",
                                 chr="all", nameChrCol="chr")


  x<-list(salmon=salmonSig$wormbaseID, dpy27=kramerDpy27$Gene_WB_ID,
          dpy21=kramerDpy21$Gene_WB_ID)
  names(x)<-c(prettyGeneName(grp), "dpy-27\n(Kramer)", "dpy-21\n(Kramer)")

  p1<-ggVennDiagram(x) + ggtitle(label=paste0(grp," vs Kramer(2015): |lfc|>", lfcVal, ", padj<",padjVal))



  lfcVal<-0.5

  kramerDpy27<-getSignificantGenes(kramer, padj=padjVal, lfc=lfcVal,
                                   namePadjCol="dpy27_RNAi_L3_padj",
                                   nameLfcCol="dpy27_RNAi_L3_log2_fold_change",
                                   direction="both",
                                   chr="all", nameChrCol="chr", outPath=".")

  kramerDpy21<-getSignificantGenes(kramer, padj=padjVal, lfc=lfcVal,
                                   namePadjCol="dpy21_mutant_L3_padj",
                                   nameLfcCol="dpy21_mutant_L3_log2_fold_change",
                                   direction="both",
                                   chr="all", nameChrCol="chr")

  salmonSig<-getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                 namePadjCol="padj",
                                 nameLfcCol="log2FoldChange",
                                 direction="both",
                                 chr="all", nameChrCol="chr")


  x<-list(salmon=salmonSig$wormbaseID, dpy27=kramerDpy27$Gene_WB_ID,
          dpy21=kramerDpy21$Gene_WB_ID)

  names(x)<-c(prettyGeneName(grp), "dpy-27\n(Kramer)", "dpy-21\n(Kramer)")

  p2<-ggVennDiagram(x) + ggtitle(label=paste0(grp," vs Kramer(2015): |lfc|>", lfcVal, ", padj<",padjVal))

  lfcVal<-0.75

  kramerDpy27<-getSignificantGenes(kramer, padj=padjVal, lfc=lfcVal,
                                   namePadjCol="dpy27_RNAi_L3_padj",
                                   nameLfcCol="dpy27_RNAi_L3_log2_fold_change",
                                   direction="both",
                                   chr="all", nameChrCol="chr", outPath=".")

  kramerDpy21<-getSignificantGenes(kramer, padj=padjVal, lfc=lfcVal,
                                   namePadjCol="dpy21_mutant_L3_padj",
                                   nameLfcCol="dpy21_mutant_L3_log2_fold_change",
                                   direction="both",
                                   chr="all", nameChrCol="chr")

  salmonSig<-getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                 namePadjCol="padj",
                                 nameLfcCol="log2FoldChange",
                                 direction="both",
                                 chr="all", nameChrCol="chr")


  x<-list(salmon=salmonSig$wormbaseID, dpy27=kramerDpy27$Gene_WB_ID,
          dpy21=kramerDpy21$Gene_WB_ID)
  names(x)<-c(prettyGeneName(grp), "dpy-27\n(Kramer)", "dpy-21\n(Kramer)")

  p3<-ggVennDiagram(x) + ggtitle(label=paste0(grp," vs Kramer(2015): |lfc|>", lfcVal, ", padj<",padjVal))


  lfcVal<-1

  kramerDpy27<-getSignificantGenes(kramer, padj=padjVal, lfc=lfcVal,
                                   namePadjCol="dpy27_RNAi_L3_padj",
                                   nameLfcCol="dpy27_RNAi_L3_log2_fold_change",
                                   direction="both",
                                   chr="all", nameChrCol="chr", outPath=".")

  kramerDpy21<-getSignificantGenes(kramer, padj=padjVal, lfc=lfcVal,
                                   namePadjCol="dpy21_mutant_L3_padj",
                                   nameLfcCol="dpy21_mutant_L3_log2_fold_change",
                                   direction="both",
                                   chr="all", nameChrCol="chr")

  salmonSig<-getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                 namePadjCol="padj",
                                 nameLfcCol="log2FoldChange",
                                 direction="both",
                                 chr="all", nameChrCol="chr")


  x<-list(salmon=salmonSig$wormbaseID, dpy27=kramerDpy27$Gene_WB_ID,
          dpy21=kramerDpy21$Gene_WB_ID)
  names(x)<-c(prettyGeneName(grp), "dpy-27\n(Kramer)", "dpy-21\n(Kramer)")

  p4<-ggVennDiagram(x) + ggtitle(label=paste0(grp," vs Kramer(2015): |lfc|>", lfcVal, ", padj<",padjVal))

  p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
  ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_bestCutoff_",grp,"_padj",
                                  formatC(padjVal,format="e",digits=0),
                                  "_lfc", "0-1",".pdf"),
                  plot=p, device="pdf",width=19,height=20,units="cm")
}





###############################
## chrX upregulated
###############################

dim(kramer)
dim(salmon)
idx<-match(kramer$Gene_WB_ID, salmon$wormbaseID)
kramer$chr<-salmon$chr[idx]
kramer<-kramer[!is.na(kramer$chr),]


padjVal=0.01

lfcVal=0
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")

idx<-!is.na(kramer$dpy27_RNAi_L3_padj) &
      kramer$dpy27_RNAi_L3_padj < padjVal &
      kramer$dpy27_RNAi_L3_log2_fold_change > lfcVal &
      kramer$chr=="chrX"
kramerdpy27dcc<-kramer[idx,]


idx<-!is.na(kramer$dpy21_mutant_L3_padj) &
  kramer$dpy21_mutant_L3_padj < padjVal &
  kramer$dpy21_mutant_L3_log2_fold_change > lfcVal &
  kramer$chr=="chrX"
kramerdpy21dcc<-kramer[idx,]

x<-list(salmon=salmondcc$wormbaseID, dpy27=kramerdpy27dcc$Gene_WB_ID,
        dpy21=kramerdpy21dcc$Gene_WB_ID)

p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))



lfcVal=0.5
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")

idx<-!is.na(kramer$dpy27_RNAi_L3_padj) &
  kramer$dpy27_RNAi_L3_padj < padjVal &
  kramer$dpy27_RNAi_L3_log2_fold_change > lfcVal &
  kramer$chr=="chrX"
kramerdpy27dcc<-kramer[idx,]


idx<-!is.na(kramer$dpy21_mutant_L3_padj) &
  kramer$dpy21_mutant_L3_padj < padjVal &
  kramer$dpy21_mutant_L3_log2_fold_change > lfcVal &
  kramer$chr=="chrX"

kramerdpy21dcc<-kramer[idx,]

x<-list(salmon=salmondcc$wormbaseID, dpy27=kramerdpy27dcc$Gene_WB_ID,
        dpy21=kramerdpy21dcc$Gene_WB_ID)

p2<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))



lfcVal=0.75
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")

idx<-!is.na(kramer$dpy27_RNAi_L3_padj) &
  kramer$dpy27_RNAi_L3_padj < padjVal &
  kramer$dpy27_RNAi_L3_log2_fold_change > lfcVal &
  kramer$chr=="chrX"
kramerdpy27dcc<-kramer[idx,]


idx<-!is.na(kramer$dpy21_mutant_L3_padj) &
  kramer$dpy21_mutant_L3_padj < padjVal &
  kramer$dpy21_mutant_L3_log2_fold_change > lfcVal &
  kramer$chr=="chrX"

kramerdpy21dcc<-kramer[idx,]

x<-list(salmon=salmondcc$wormbaseID, dpy27=kramerdpy27dcc$Gene_WB_ID,
        dpy21=kramerdpy21dcc$Gene_WB_ID)

p3<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))



lfcVal=1
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")

idx<-!is.na(kramer$dpy27_RNAi_L3_padj) &
  kramer$dpy27_RNAi_L3_padj < padjVal &
  kramer$dpy27_RNAi_L3_log2_fold_change > lfcVal &
  kramer$chr=="chrX"
kramerdpy27dcc<-kramer[idx,]


idx<-!is.na(kramer$dpy21_mutant_L3_padj) &
  kramer$dpy21_mutant_L3_padj < padjVal &
  kramer$dpy21_mutant_L3_log2_fold_change > lfcVal &
  kramer$chr=="chrX"

kramerdpy21dcc<-kramer[idx,]

x<-list(salmon=salmondcc$wormbaseID, dpy27=kramerdpy27dcc$Gene_WB_ID,
        dpy21=kramerdpy21dcc$Gene_WB_ID)

p4<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_XchrUp_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-1",".pdf"),
                plot=p, device="pdf",width=19,height=20,units="cm")



###############################
## classical DCC genes
###############################

pubDCC<-readRDS("/Users/semple/Documents/MeisterLab/dSMF/DCgenes/published_DCCgr.rds")
pubNDC<-readRDS("/Users/semple/Documents/MeisterLab/dSMF/DCgenes/published_NDCgr.rds")


padjVal=0.01
lfcVal=0
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, DC=pubDCC$wormbaseID,
        nonDC=pubNDC$wormbaseID)
p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.5
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, DC=pubDCC$wormbaseID,
        nonDC=pubNDC$wormbaseID)
p2<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.75
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, DC=pubDCC$wormbaseID,
        nonDC=pubNDC$wormbaseID)
p3<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))

lfcVal=1
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, DC=pubDCC$wormbaseID,
        nonDC=pubNDC$wormbaseID)
p4<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_papers_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-1",".pdf"),
                plot=p, device="pdf",width=19,height=20,units="cm")


padjVal=0.05
lfcVal=0
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, DC=pubDCC$wormbaseID,
        nonDC=pubNDC$wormbaseID)
p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.5
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, DC=pubDCC$wormbaseID,
        nonDC=pubNDC$wormbaseID)
p2<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.75
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, DC=pubDCC$wormbaseID,
        nonDC=pubNDC$wormbaseID)
p3<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))

lfcVal=1
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, DC=pubDCC$wormbaseID,
        nonDC=pubNDC$wormbaseID)
p4<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_papers_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-1",".pdf"),
                plot=p, device="pdf",width=19,height=20,units="cm")




###############################
## Jans 2009
###############################

JansDC<-readRDS("/Users/semple/Documents/MeisterLab/dSMF/DCgenes/Jans2009_DCgr.rds")
JansNDC<-readRDS("/Users/semple/Documents/MeisterLab/dSMF/DCgenes/Jans2009_NDCgr.rds")


padjVal=0.01
lfcVal=0
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, JansDC=JansDC$wormbaseID,
        JansNDC=JansNDC$wormbaseID)
p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.5
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, JansDC=JansDC$wormbaseID,
        JansNDC=JansNDC$wormbaseID)
p2<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.75
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, JansDC=JansDC$wormbaseID,
        JansNDC=JansNDC$wormbaseID)
p3<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))

lfcVal=1
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, JansDC=JansDC$wormbaseID,
        JansNDC=JansNDC$wormbaseID)
p4<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_Jans2009_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-1",".pdf"),
                plot=p, device="pdf",width=19,height=20,units="cm")


padjVal=0.05
lfcVal=0
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, JansDC=JansDC$wormbaseID,
        JansNDC=JansNDC$wormbaseID)
p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.5
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, JansDC=JansDC$wormbaseID,
        JansNDC=JansNDC$wormbaseID)
p2<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.75
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, JansDC=JansDC$wormbaseID,
        JansNDC=JansNDC$wormbaseID)
p3<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))

lfcVal=1
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbaseID, JansDC=JansDC$wormbaseID,
        JansNDC=JansNDC$wormbaseID)
p4<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_Jans2009_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-1",".pdf"),
                plot=p, device="pdf",width=19,height=20,units="cm")


###############################
## Jans 2009 vs Kramer
###############################

kramer<-read_excel(kramerFileName,col_types=c(rep("text",3),rep("numeric",30)))
names(kramer)
kramer<-kramer[,c(1:3,grep("_L3_",names(kramer)))]

padjVal<-0.05
lfcVal<-0
idx<-!is.na(kramer$dpy27_RNAi_L3_padj) & abs(kramer$dpy27_RNAi_L3_log2_fold_change)>lfcVal & kramer$dpy27_RNAi_L3_padj<padjVal
kramerDpy27<-kramer[idx,]

idx<-!is.na(kramer$dpy21_mutant_L3_padj) & abs(kramer$dpy21_mutant_L3_log2_fold_change)>lfcVal & kramer$dpy21_mutant_L3_padj<padjVal
kramerDpy21<-kramer[idx,]


x<-list(JansDC=JansDC$wormbaseID, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p1<-ggVennDiagram(x) + ggtitle(label=paste0("Jans DC vs Kramer(2015): lfc>", lfcVal, ", padj<",padjVal))


x<-list(JansNDC=JansNDC$wormbaseID, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p2<-ggVennDiagram(x) + ggtitle(label=paste0("Jans NDC vs Kramer(2015): lfc>", lfcVal, ", padj<",padjVal))


lfcVal<-0.5
idx<-!is.na(kramer$dpy27_RNAi_L3_padj) & abs(kramer$dpy27_RNAi_L3_log2_fold_change)>lfcVal & kramer$dpy27_RNAi_L3_padj<padjVal
kramerDpy27<-kramer[idx,]

idx<-!is.na(kramer$dpy21_mutant_L3_padj) & abs(kramer$dpy21_mutant_L3_log2_fold_change)>lfcVal & kramer$dpy21_mutant_L3_padj<padjVal
kramerDpy21<-kramer[idx,]


x<-list(JansDC=JansDC$wormbaseID, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p3<-ggVennDiagram(x) + ggtitle(label=paste0("Jans DC vs Kramer(2015): lfc>", lfcVal, ", padj<",padjVal))



x<-list(JansNDC=JansNDC$wormbaseID, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p4<-ggVennDiagram(x) + ggtitle(label=paste0("Jans NDC vs Kramer(2015): lfc>", lfcVal, ", padj<",padjVal))


p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_Jans2009vKramer2015_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-0.5",".pdf"),
                plot=p, device="pdf",width=19,height=20,units="cm")


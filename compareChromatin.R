library(rtracklayer)
library(ggplot2)
library(BSgenome.Celegans.UCSC.ce11)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

source("functions.R")
source("./variableSettings.R")

scriptName <- "compareChromatin"
print(scriptName)

if(filterData){
  fileNamePrefix<-filterPrefix
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
} else {
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
}

makeDirs(outPath,dirNameList=paste0(c("plots/","tracks/"),
                                    paste0(dirname(fileNamePrefix),"/",
                                           scriptName)))

############################-
# overlap of states with significant genes---------
############################-

states<-import.bed("./publicData/chromStates_L3_Evans2016_ce11.bed")
#domains<-import.bed("./publicData/chromDomains_L3_Evans2016_ce11.bed")
#seqlevels(domains)<-seqlevels(states)

stateClrs<-c("#fe0003","#f59745","#008100","#74943a",
             "#c4d69c","#05ff7f","#ceff65","#fd0082",
             "#ff70d0","#ffb6c4","#7f00ff","#1900fe",
             "#528dd4","#814006","#b8cde4","#808080",
             "#dcdac3","#c4bc97","#938b54","#141414")

getStateOLtable<-function(gr,states){
  ol<-data.frame(findOverlaps(gr,states,ignore.strand=T,minoverlap=10))
  gr$name<-paste0(seqnames(gr),":",start(gr),"-",end(gr))
  ol$name<-gr$name[ol$queryHits]
  ol$XvA<-ifelse(unlist(strsplit(ol$name,":.*$"))=="chrX","X","A")
  ol$state<-factor(states$score[ol$subjectHits],levels=1:20)
  ol$width<-width(states)[ol$subjectHits]
  df<-ol%>%group_by(state,XvA) %>% summarise(stateFrequency=n(),stateWidth=sum(width))
  allStates<-data.frame(state=factor(1:20, levels=1:20))
  df<-left_join(allStates,df)
  df[is.na(df)]<-0
  return(df)
}


## significant genes-----
sigTables<-list()
plotList<-list()
grp<-useContrasts[1]
for (grp in useContrasts[c(1,3,4,6,7,8,9)]){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  up<-GRanges(getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                    namePadjCol="padj",
                                    nameLfcCol="log2FoldChange",
                                    direction="gt",
                                    chr="all", nameChrCol="chr"))
  stateTblup<-getStateOLtable(up,states)
  stateTblup$upVdown<-"up"
  stateTblup$SMC<-grp

  down<-GRanges(getSignificantGenes(salmon, padj=padjVal, lfc=-lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="lt",
                        chr="all", nameChrCol="chr"))
  stateTbldown<-getStateOLtable(down,states)
  stateTbldown$upVdown<-"down"
  stateTbldown$SMC<-grp
  sigTables[[grp]]<-rbind(stateTblup,stateTbldown)
}


sig<-do.call(rbind,sigTables)
row.names(sig)<-NULL
sig$SMC<-factor(sig$SMC,levels=useContrasts[c(1,3,4,6,7,8,9)])

p1<-ggplot(sig[sig$XvA=="A",],aes(x=state,y=stateFrequency,fill=state))+
  geom_bar(stat="identity") + facet_grid(rows=vars(SMC),cols=vars(upVdown))+
  theme_bw() + theme(legend.position="none") +
  ggtitle(label="Counts of states overlapping significant genes on autosomes") +
  scale_fill_manual(values=stateClrs)
p1

p2<-ggplot(sig[sig$XvA=="A",],aes(x=state,y=stateWidth,fill=state))+
  geom_bar(stat="identity") + facet_grid(rows=vars(SMC),cols=vars(upVdown))+
  theme_bw() + theme(legend.position="none") +
  ggtitle(label="Width of states overlapping significant genes on autosomes") +
  scale_fill_manual(values=stateClrs)
p2


p3<-ggplot(sig[sig$XvA=="X",],aes(x=state,y=stateFrequency,fill=state))+
  geom_bar(stat="identity") + facet_grid(rows=vars(SMC),cols=vars(upVdown))+
  theme_bw() + theme(legend.position="none") +
  ggtitle(label="Counts of states overlapping significant genes on chrX") +
  scale_fill_manual(values=stateClrs)
p3

p4<-ggplot(sig[sig$XvA=="X",],aes(x=state,y=stateWidth,fill=state))+
  geom_bar(stat="identity") + facet_grid(rows=vars(SMC),cols=vars(upVdown))+
  theme_bw() + theme(legend.position="none") +
  ggtitle(label="Width of states overlapping significant genes on chrX") +
  scale_fill_manual(values=stateClrs)
p4


p<-ggarrange(p1,p2,p3,p4,ncol=1,nrow=1)
ggexport(p,filename=paste0(outPath,"/plots/",outputNamePrefix,"chromStates_p",padjVal,"_lfc",lfcVal,".pdf"),height=9,width=9,units="cm",device="pdf")


## any lfc genes-----
sigTables<-list()
plotList<-list()
localLFCval<-0.3
grp<-useContrasts[1]
for (grp in useContrasts[c(1,3,4,6,7,8,9)]){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  up<-GRanges(getSignificantGenes(salmon, padj=padjVal, lfc=localLFCval,
                                  namePadjCol="padj",
                                  nameLfcCol="log2FoldChange",
                                  direction="gt",
                                  chr="all", nameChrCol="chr"))
  stateTblup<-getStateOLtable(up,states)
  stateTblup$upVdown<-"up"
  stateTblup$SMC<-grp

  down<-GRanges(getSignificantGenes(salmon, padj=padjVal, lfc=-localLFCval,
                                    namePadjCol="padj",
                                    nameLfcCol="log2FoldChange",
                                    direction="lt",
                                    chr="all", nameChrCol="chr"))
  stateTbldown<-getStateOLtable(down,states)
  stateTbldown$upVdown<-"down"
  stateTbldown$SMC<-grp
  sigTables[[grp]]<-rbind(stateTblup,stateTbldown)
}


sig<-do.call(rbind,sigTables)
row.names(sig)<-NULL
sig$SMC<-factor(sig$SMC,levels=useContrasts[c(1,3,4,6,7,8,9)])

p1<-ggplot(sig[sig$XvA=="A",],aes(x=state,y=stateFrequency,fill=state))+
  geom_bar(stat="identity") + facet_grid(rows=vars(SMC),cols=vars(upVdown))+
  theme_bw() + theme(legend.position="none") +
  ggtitle(label="Counts of states overlapping significant genes on autosomes") +
  scale_fill_manual(values=stateClrs)
p1

p2<-ggplot(sig[sig$XvA=="A",],aes(x=state,y=stateWidth,fill=state))+
  geom_bar(stat="identity") + facet_grid(rows=vars(SMC),cols=vars(upVdown))+
  theme_bw() + theme(legend.position="none") +
  ggtitle(label="Width of states overlapping significant genes on autosomes") +
  scale_fill_manual(values=stateClrs)
p2


p3<-ggplot(sig[sig$XvA=="X",],aes(x=state,y=stateFrequency,fill=state))+
  geom_bar(stat="identity") + facet_grid(rows=vars(SMC),cols=vars(upVdown))+
  theme_bw() + theme(legend.position="none") +
  ggtitle(label="Counts of states overlapping significant genes on chrX") +
  scale_fill_manual(values=stateClrs)
p3

p4<-ggplot(sig[sig$XvA=="X",],aes(x=state,y=stateWidth,fill=state))+
  geom_bar(stat="identity") + facet_grid(rows=vars(SMC),cols=vars(upVdown))+
  theme_bw() + theme(legend.position="none") +
  ggtitle(label="Width of states overlapping significant genes on chrX") +
  scale_fill_manual(values=stateClrs)
p4


p<-ggarrange(p1,p2,p3,p4,ncol=1,nrow=1)
ggexport(p,filename=paste0(outPath,"/plots/",outputNamePrefix,"chromStates_p",padjVal,"_lfc",localLFCval,".pdf"),height=9,width=9,units="cm",device="pdf")


### all genes
sigTables<-list()
plotList<-list()
grp<-useContrasts[1]
for (grp in useContrasts[c(1,3,4,6,7,8,9)]){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  gr<-GRanges(salmon[!is.na(salmon$padj),])
  stateTblgr<-getStateOLtable(gr,states)
  stateTblgr$type<-"tested"
  stateTblgr$SMC<-grp

  gr1<-GRanges(salmon)
  stateTblgr1<-getStateOLtable(gr1,states)
  stateTblgr1$type<-"all"
  stateTblgr1$SMC<-grp
  sigTables[[grp]]<-rbind(stateTblgr,stateTblgr1)
}


sig<-do.call(rbind,sigTables)
row.names(sig)<-NULL
sig$SMC<-factor(sig$SMC,levels=useContrasts[c(1,3,4,6,7,8,9)])

p1<-ggplot(sig[sig$XvA=="A",],aes(x=state,y=stateFrequency,fill=state))+
  geom_bar(stat="identity") + facet_grid(rows=vars(SMC),cols=vars(type))+
  theme_bw() + theme(legend.position="none") +
  ggtitle(label="Counts of states overlapping all genes on autosomes") +
  scale_fill_manual(values=stateClrs)
p1

p2<-ggplot(sig[sig$XvA=="A",],aes(x=state,y=stateWidth,fill=state))+
  geom_bar(stat="identity") + facet_grid(rows=vars(SMC),cols=vars(type))+
  theme_bw() + theme(legend.position="none") +
  ggtitle(label="Width of states overlapping all genes on autosomes") +
  scale_fill_manual(values=stateClrs)
p2


p3<-ggplot(sig[sig$XvA=="X",],aes(x=state,y=stateFrequency,fill=state))+
  geom_bar(stat="identity") + facet_grid(rows=vars(SMC),cols=vars(type))+
  theme_bw() + theme(legend.position="none") +
  ggtitle(label="Counts of states overlapping all genes on chrX") +
  scale_fill_manual(values=stateClrs)
p3

p4<-ggplot(sig[sig$XvA=="X",],aes(x=state,y=stateWidth,fill=state))+
  geom_bar(stat="identity") + facet_grid(rows=vars(SMC),cols=vars(type))+
  theme_bw() + theme(legend.position="none") +
  ggtitle(label="Width of states overlapping all genes on chrX") +
  scale_fill_manual(values=stateClrs)
p4


p<-ggarrange(p1,p2,p3,p4,ncol=1,nrow=1)
ggexport(p,filename=paste0(outPath,"/plots/",outputNamePrefix,"chromStates_all.pdf"),height=9,width=9,units="cm",device="pdf")

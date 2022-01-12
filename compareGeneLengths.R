
library(lattice)
library(ggplot2)
library(dplyr)

source("functions.R")
source("./variableSettings.R")
scriptName <- "compareGeneLengths"
print(scriptName)

if(filterData){
  fileNamePrefix<-filterPrefix
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
} else {
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
}

makeDirs(outPath,dirNameList=paste0(c("plots/"),
                                    paste0(dirname(fileNamePrefix),"/",
                                           scriptName)))


###########################
## compare gene lengths
##########################

getSimSig<-function(simTbl,numSims,grp, chrSubset,
                    metricNames=c("min","mean","median","max")){
  recordSig<-list()
  for(metricName in metricNames){
    qval=sum(simTbl[numSims+1,metricName]>simTbl[1:numSims,metricName])/numSims
    recordSig[[paste(grp,chrSubset,metricName,sep="_")]]<-data.frame(group=grp,
                                                                     subset=chrSubset,
                                                                     metric=metricName,
                                                                     quantile=qval)
  }
  return(do.call(rbind,recordSig))
}

geneLengthHist<-function(simTbl,metricName,numSims,grp,chrSubset,direction){
  qval=sum(simTbl[numSims+1,metricName]>simTbl[1:numSims,metricName])/numSims
  hist(simTbl[1:numSims,metricName], main=paste0(grp," ",chrSubset,"(",direction,
                                                 "): Sim ",metricName," gene length"),
       sub=paste0("quantile:",qval,ifelse(qval<0.05 | qval>0.95,"**","")),
       xlab=paste0(metricName," gene length"), xlim=c(min(simTbl[,metricName])*0.9,
                                                      max(simTbl[,metricName])*1.1),
       breaks=50)
  abline(v=simTbl[numSims+1,metricName],col="red")
  abline(v=median(simTbl[1:numSims,metricName]),lty=2,col="grey40")
  legend("topright",legend=c("true value","sims median"),lty=c(1,2),col=c("red","black"))
}

doSims<-function(numSims,sig.gr,bg.gr,grp,chrSubset,padj,lfc,direction,outputNamePrefix,outPath="."){
  #get a randomly sampled distribution of same size
  simTbl<-as.data.frame(matrix(NA,nrow=numSims+1,ncol=5))
  #names(simTbl)<-c("min","mean","median","max","meanT10")
  names(simTbl)<-c("mean","median")
  for(i in 1:numSims){
    sampleIdxs<-sample(1:length(bg.gr),length(sig.gr))
    sampledGenes<-bg.gr[sampleIdxs]
    geneLengths<-width(sampledGenes[order(width(sampledGenes),decreasing=T)])
    #simTbl[i,"min"]<-min(geneLengths)
    simTbl[i,"mean"]<-mean(geneLengths)
    simTbl[i,"median"]<-median(geneLengths)
    #simTbl[i,"max"]<-max(geneLengths)
    #simTbl[i,"meanT10"]<-mean(geneLengths[1:10])
  }
  geneLengths<-width(sig.gr[order(width(sig.gr),decreasing=T)])
  #simTbl[numSims+1,"min"]<-min(geneLengths)
  simTbl[numSims+1,"mean"]<-mean(geneLengths)
  simTbl[numSims+1,"median"]<-median(geneLengths)
  #simTbl[numSims+1,"max"]<-max(geneLengths)
  #simTbl[numSims+1,"meanT10"]<-mean(geneLengths[1:10])
  pdf(file=paste0(outPath, "/plots/",outputNamePrefix,"geneLengthHist_",
                  grp,"_",chrSubset,"_",direction,"_padj",padj,"_lfc", lfc,".pdf"),
      width=11, height=3.5, paper="a4r")
  par(mfrow=c(1,3))
  hist(geneLengths, main=paste0(grp," ",chrSubset,": Gene length distribution"),
       sub=paste0(length(sig.gr)," regulated genes, ",length(bg.gr)," ",direction," genes"),
       xlab="Gene length",col="lightblue",breaks=50)
  #geneLengthHist(simTbl=simTbl,metricName="min",numSims=numSims,grp=grp,chrSubset=chrSubset)
  geneLengthHist(simTbl=simTbl,metricName="mean",numSims=numSims,grp=grp,chrSubset=chrSubset,direction=direction)
  geneLengthHist(simTbl=simTbl,metricName="median",numSims=numSims,grp=grp,chrSubset=chrSubset,direction=direction)
  #geneLengthHist(simTbl=simTbl,metricName="max",numSims=numSims,grp=grp,chrSubset=chrSubset)
  #geneLengthHist(simTbl=simTbl,metricName="meanT10",numSims=numSims,grp=grp,chrSubset=chrSubset)
  dev.off()
  recordSig<-getSimSig(simTbl,numSims,grp, chrSubset,metricNames=c("mean","median"))
  #recordSig<-getSimSig(simTbl,numSims,grp, chrSubset,metricNames=c("min","mean","median","max"))
  return(recordSig)
}


# compare lengths of genes
sigTable<-list()
bgTable<-list()
sigGR<-list()
bgGR<-list()
numSims<-1000
recordSig<-list()
set.seed(97232151)
#df<-data.frame(matrix(nrow=10,ncol=3))
#names(df)<-useContrasts

compareLengthToSims<-function(grp,contrastNames,fileNamePrefix,padjVal,lfcVal,direction) {
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  bgTable<-data.frame(getSignificantGenes(salmon, padj=1, lfc=0,
                               namePadjCol="padj",
                               nameLfcCol="log2FoldChange",
                               direction=direction,
                               chr="all", nameChrCol="chr"))
  sigTable<-data.frame(getSignificantGenes(salmon, padj=padjVal,
                                lfc=ifelse(direction=="lt",-lfcVal,lfcVal),
                                namePadjCol="padj",
                                nameLfcCol="log2FoldChange",
                                direction=direction,
                                chr="all", nameChrCol="chr"))
  sigTable<-sigTable[!is.na(sigTable$chr),] # removes mtDNA genes
  bgTable<-bgTable[!is.na(bgTable$chr),] # removes mtDNA genes
  sigGR<-GRanges(seqnames=sigTable$chr,
                 ranges=IRanges(start=sigTable$start,
                                end=sigTable$end),
                 strand=sigTable$strand)
  bgGR<-GRanges(seqnames=bgTable$chr,
                ranges=IRanges(start=bgTable$start,
                               end=bgTable$end),
                strand=bgTable$strand)
  mcols(sigGR)<-sigTable
  mcols(bgGR)<-bgTable


  # check if datasets have chrX genes included
  includeChrX<-"chrX" %in% sigTable$chr
  if(includeChrX){
    chrXgr<-sigGR[sigGR$chr=="chrX"]
    chrXgrBg<-bgGR[bgGR$chr=="chrX"]
    if(length(chrXgr)>20){
      recordSig[[paste(grp,"chrX",sep="_")]]<-doSims(numSims=numSims,sig.gr=chrXgr,
                                                     bg.gr=chrXgrBg,grp=grp,
                                                     chrSubset="chrX",
                                                     direction=direction,
                                                     padj=padjVal,lfc=lfcVal,
                                                     outputNamePrefix=outputNamePrefix,
                                                     outPath=".")
    }
  }

  chrAgr<-sigGR[sigGR$chr!="chrX"]
  chrAgrBg<-bgGR[bgGR$chr!="chrX"]
  if(length(chrAgr)>20){
    recordSig[[paste(grp,"chrA",sep="_")]]<-doSims(numSims=numSims,sig.gr=chrAgr,
                                                   bg.gr=chrAgrBg,grp=grp,
                                                   chrSubset="chrA",
                                                   direction=direction,
                                                   padj=padjVal,lfc=lfcVal,
                                                   outputNamePrefix=outputNamePrefix,
                                                   outPath=".")
  }
  return(recordSig)
}


for (grp in useContrasts){
  print(grp)
  print("both")
  compareLengthToSims(grp,contrastNames,fileNamePrefix,padjVal,lfcVal,direction="both")
  print("gt")
  compareLengthToSims(grp,contrastNames,fileNamePrefix,padjVal,lfcVal,direction="gt")
  print("lt")
  compareLengthToSims(grp,contrastNames,fileNamePrefix,padjVal,lfcVal,direction="lt")
}
# par(mfrow=c(1,1))
# recordSig1<-do.call(rbind,recordSig)
# rownames(recordSig1)<-NULL
# recordSig1 %>%  filter(metric=="mean") %>% group_by(subset)
# recordSig1 %>%  filter(metric=="median") %>% group_by(subset)





############# plot length vs basal expression
chrSubset="autosomes"
localPadj=0.05
localLFC=0
grp=useContrasts[1]
listTbls<-list()
for(grp in useContrasts[c(3,6,7,8,9)]){
  salmon<-data.frame(readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds")))
  sig<-getSignificantGenes(salmon, padj=localPadj, lfc=localLFC,
                                        namePadjCol="padj",
                                        nameLfcCol="log2FoldChange",
                                        direction="both",
                                        chr=chrSubset, nameChrCol="chr")
  sig$geneLength<-sig$end-sig$start
  sig$upVdown<-factor(ifelse(sig$log2FoldChange>0,"up","down"))
  sig$SMC<-grp
  listTbls[[grp]]<-sig
}


allSig<-do.call(rbind,listTbls)
allSig$SMC<-factor(allSig$SMC,levels=useContrasts[c(1,3,4,6,7,8,9)])
p1<-ggplot(allSig,aes(x=log2(geneLength),y=log2(baseMean),color=log2FoldChange)) +
    geom_point(size=0.4) +
    scale_color_gradient2(low=muted("#ff000055"),mid="#ffffff22",
                          high=muted("#0000ff55"), na.value="#ffffff22",
                          limits=c(-0.5,0.5),oob=scales::squish,name="Log2FC")+
    facet_grid(rows=vars(SMC)) +theme_bw()+
  ggtitle(paste0("Significantly changed genes on ",chrSubset," LFC>",localLFC))+
  theme(legend.position = "bottom")



uniqGenes<-allSig %>% distinct(wormbaseID,geneLength)

allSig$lengthBin<-cut(allSig$geneLength,quantile(uniqGenes$geneLength,seq(0,1,0.1)),
                      dig.lab=0,ordered_result=T,right=T,include.lowest=T)


labs<-data.frame(lower = factor( as.numeric(sub("(\\(|\\[)(.+),.*", "\\2", levels(allSig$lengthBin) ))),
                 upper = factor( as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", levels(allSig$lengthBin) ))))

levels(allSig$lengthBin)<-paste(levels(labs$lower), levels(labs$upper),sep="-")
#allSig[is.na(allSig$lengthBin),]
p2<-ggplot(allSig, aes(x=lengthBin,fill=upVdown)) + geom_bar(position="fill") +
  facet_grid(rows=vars(SMC)) + theme_bw() +
  scale_fill_manual(values=c(muted("red"),muted("blue")),name="Log2FC") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "bottom") +
  ylab("Fraction of genes") +
  ggtitle(paste0("Proportion up vs down by gene length"))
#+
#  geom_text(aes(label = ..count.., y=..count..),  stat = "count",
#            position = position_fill(vjust = .5),col="white",angle=90)

#p2

p<-ggarrange(p1,p2,ncol=2,widths=c(3,1.5))
ggsave(filename=paste0(outPath, "/plots/",outputNamePrefix,"geneLengthBaseMean&LFC_",
                       chrSubset,"_padj",localPadj,"_lfc", localLFC,".pdf"), plot=p,
       width=8, height=11,device="pdf")



### binned expression
allSig$expressionBin<-cut(allSig$baseMean,quantile(allSig$baseMean,seq(0,1,0.1)),
                      dig.lab=0,ordered_result=T,right=T,include.lowest=T)

allSig %>% dplyr::group_by(SMC,expressionBin,upVdown) %>% dplyr::summarise(count=n())

labs<-data.frame(lower = factor( as.numeric(sub("(\\(|\\[)(.+),.*", "\\2", levels(allSig$expressionBin) ))),
                 upper = factor( as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", levels(allSig$expressionBin) ))))

levels(allSig$expressionBin)<-paste(levels(labs$lower), levels(labs$upper),sep="-")
#allSig[is.na(allSig$expressionBin),]
p3<-ggplot(allSig, aes(x=expressionBin,fill=upVdown)) + geom_bar(position="fill") +
  facet_grid(rows=vars(SMC)) + theme_bw() +
  scale_fill_manual(values=c(muted("red"),muted("blue")),name="Log2FC") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "bottom") +
  ylab("Fraction of genes") +
  ggtitle(paste0("Proportion up vs down by gene length"))
p3

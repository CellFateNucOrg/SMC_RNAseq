
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
                                    paste0("p",padjVal,"_lfc",lfcVal,"/",
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

geneLengthHist<-function(simTbl,metricName,numSims,grp,chrSubset){
  qval=sum(simTbl[numSims+1,metricName]>simTbl[1:numSims,metricName])/numSims
  hist(simTbl[1:numSims,metricName], main=paste0(grp," ",chrSubset,": Sim ",metricName," gene length"),
       sub=paste0("quantile:",qval,ifelse(qval<0.05 | qval>0.95,"**","")),
       xlab=paste0(metricName," gene length"), xlim=c(min(simTbl[,metricName])*0.9,
                                                      max(simTbl[,metricName])*1.1),
       breaks=50)
  abline(v=simTbl[numSims+1,metricName],col="red")
  abline(v=median(simTbl[1:numSims,metricName]),lty=2,col="grey40")
  legend("topright",legend=c("true value","sims median"),lty=c(1,2),col=c("red","black"))
}

doSims<-function(numSims,sig.gr,bg.gr,grp,chrSubset,padj,lfc,outputNamePrefix,outPath="."){
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
                  grp,"_",chrSubset,"_padj",padj,"_lfc", lfc,".pdf"),
      width=11, height=3.5, paper="a4r")
  par(mfrow=c(1,3))
  hist(geneLengths, main=paste0(grp," ",chrSubset,": Gene length distribution"),
       sub=paste0(length(sig.gr)," regulated genes, ",length(bg.gr)," all expressed genes"),
       xlab="Gene length",col="lightblue",breaks=50)
  #geneLengthHist(simTbl=simTbl,metricName="min",numSims=numSims,grp=grp,chrSubset=chrSubset)
  geneLengthHist(simTbl=simTbl,metricName="mean",numSims=numSims,grp=grp,chrSubset=chrSubset)
  geneLengthHist(simTbl=simTbl,metricName="median",numSims=numSims,grp=grp,chrSubset=chrSubset)
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
for (grp in names(contrastNames)){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

  bgTable<-as.data.frame(salmon[!is.na(salmon$padj),])
  sigTable<-as.data.frame(getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                              namePadjCol="padj",
                                              nameLfcCol="log2FoldChange",
                                              direction="both",
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
                                                   padj=padjVal,lfc=lfcVal,
                                                   outputNamePrefix=outputNamePrefix,
                                                   outPath=".")
  }

}
par(mfrow=c(1,1))
recordSig1<-do.call(rbind,recordSig)
rownames(recordSig1)<-NULL
recordSig1 %>%  filter(metric=="mean") %>% group_by(subset)
recordSig1 %>%  filter(metric=="median") %>% group_by(subset)

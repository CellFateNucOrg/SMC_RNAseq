library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(BSgenome.Celegans.UCSC.ce11)
library(RColorBrewer)
library(seqplots)
library(tidyr)

source("functions.R")
source("./variableSettings.R")

scriptName <- "compareToRexMex"
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

##################-
## distance of rex/mex from significant genes
##################-

sigGR<-list()
for (grp in useContrasts){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],
                         "_DESeq2_fullResults_p",padjVal,".rds"))
  sigGR[[grp]]<-GRanges(getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                                      namePadjCol="padj",
                                                      nameLfcCol="log2FoldChange",
                                                      direction="both",
                                                      chr="all", nameChrCol="chr"))
}

dpy26bw<-sigGR[["dpy26"]]
seqinfo(dpy26bw)<-seqinfo(Celegans)
mcols(dpy26bw)$score<-dpy26bw$log2FoldChange
sdc3bw<-sigGR[["aux_sdc3BG"]]
seqinfo(sdc3bw)<-seqinfo(Celegans)
mcols(sdc3bw)$score<-sdc3bw$log2FoldChange
export.bed(dpy26bw,"dpy26sig.bed")
export.bed(sdc3bw,"sdc3sig.bed")


pdf(file=paste0(outPath, "/plots/",outputNamePrefix,"distance2rex_",
                "_padj",padjVal, "_lfc", lfcVal,".pdf"),
    width=11,height=8,paper="a4r")

# rex Motifs
rexMotifs<-import.bed(paste0(outPath,"/publicData/rexMotifs_Albritton2017_ce11.bed"))

for(grp in useContrasts){
  dist2rex<-distanceToNearest(sigGR[[grp]],rexMotifs,ignore.strand=T)
  sigGR[[grp]]$dist2rex<-NA
  sigGR[[grp]]$dist2rex[queryHits(dist2rex)]<-mcols(dist2rex)$distance
  sigGR[[grp]]$SMC<-grp
}

dfall<-do.call(rbind,lapply(sigGR,data.frame))
row.names(dfall)<-NULL
dfall$XvA<-ifelse(dfall$seqnames=="chrX","X","A")
dfall$upVdown<-ifelse(dfall$log2FoldChange>0,"up","down")
head(dfall)

p1<-ggplot(dfall,aes(x=dist2rex/1000,fill=XvA))+geom_histogram() +facet_grid(XvA~SMC) +
  theme_bw() +xlab("Distance to rex Motif (kb)") +
  ggtitle("Distance of significantly changed genes to nearest Rex Motif")
print(p1)

p2<-ggplot(dfall[dfall$XvA=="X",],aes(x=dist2rex/1000,fill=upVdown))+geom_histogram() +facet_grid(upVdown~SMC) +
  theme_bw() +xlab("Distance to rex Motif (kb)") +
  ggtitle("Distance of chrX significantly changed genes to nearest Rex Motif")
print(p2)


#### rex sites
rexSites<-import.bed(paste0(outPath,"/publicData/rexSites_Albritton2017_ce11.bed"))

for(grp in useContrasts){
  dist2rex<-distanceToNearest(sigGR[[grp]],rexSites,ignore.strand=T)
  sigGR[[grp]]$dist2rex<-NA
  sigGR[[grp]]$dist2rex[queryHits(dist2rex)]<-mcols(dist2rex)$distance
  sigGR[[grp]]$SMC<-grp
}

dfall<-do.call(rbind,lapply(sigGR,data.frame))
row.names(dfall)<-NULL
dfall$XvA<-ifelse(dfall$seqnames=="chrX","X","A")
dfall$upVdown<-ifelse(dfall$log2FoldChange>0,"up","down")
head(dfall)

# p3<-ggplot(dfall,aes(x=dist2rex/1000,fill=XvA))+geom_histogram() +facet_grid(XvA~SMC) +
#   theme_bw() +xlab("Distance to rex Site (kb)")
# print(p3)

p4<-ggplot(dfall[dfall$XvA=="X",],aes(x=dist2rex/1000,fill=upVdown))+geom_histogram() +facet_grid(upVdown~SMC) +
  theme_bw() +xlab("Distance to rex Site (kb)") +
  ggtitle("Distance of chrX significantly changed genes to nearest Rex Site")
print(p4)

dev.off()




###############-
## seqplots for rex/mex sites
###############-

# anchors - seqplots heatmaps ---------------------------------------------
####
## anchors
####

# rex Motifs
rexMotifs<-import.bed(paste0(outPath,"/publicData/rexMotifs_Albritton2017_ce11.bed"))

rexMotifs<-rexMotifs[seqnames(rexMotifs)=="chrX"]

rexMotifs<-resize(rexMotifs,width=10000,fix="center")

flankSize<-60000

smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
                  useContrasts,"_lfc.bw")
if(plotPDFs==T){
  pdf(file=paste0(outPath,"/plots/",outputNamePrefix,"rexMotifs-chrX_flank",
                  flankSize/1000,"kb.pdf"), width=19,
      height=16, paper="a4")
}

if(plotPDFs==F){
  png(filename=paste0(outPath,"/plots/",outputNamePrefix,"rexMotifs-chrX_flank",
                      flankSize/1000,"kb.png"),width=19,
      height=16, units="cm", res=150)
}

#############################
##### subsitute for getREF function from seqplots thaht has an unfixed bug.
##### Fix comes form  https://github.com/Przemol/seqplots/issues/58
#' Get reference genome
#'
#' @param genome The filename of FASTA file or genome code for BSgenome
#'
#' @return \code{DNAStringSet}
#'
#' @export
#'
getREF <- function(genome) {
  if( file.exists(file.path(Sys.getenv('root'), 'genomes', genome)) ) {
    REF <- Biostrings::readDNAStringSet( file.path(Sys.getenv('root'), 'genomes', genome) )
    names(REF) <- gsub(' .+', '', names(REF))
  } else {

    GENOMES <- BSgenome::installed.genomes(
      splitNameParts=TRUE)$genome
    if( length(GENOMES) )
      names(GENOMES) <- gsub('^BSgenome.', '', BSgenome::installed.genomes())
    if( !length(GENOMES) ) stop('No genomes installed!')

    pkg <- paste0('BSgenome.', names(GENOMES[GENOMES %in% genome]))[[1]]
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE, quietly=TRUE)
    )
    REF <- get(pkg)
  }
  return(REF)
}

assignInNamespace("getREF",getREF,ns="seqplots")
#####################


p<-getPlotSetArray(tracks=c(smcRNAseq),
                   features=c(rexMotifs),
                   refgenome="ce11", bin=100L, xmin=flankSize,
                   xmax=flankSize, type="af",
                   xanchored=10000)


dd<-plotHeatmap(p,plotz=F)
heatmapQuantiles<-sapply(dd$HLST,quantile,c(0.05,0.95),na.rm=T)
roworder<-rev(order(lapply(dd$HLST,rowSums,na.rm=T)$X1))
minVal<-min(heatmapQuantiles[1,])
maxVal<-max(heatmapQuantiles[2,])
plotHeatmap(p,main="chrX rex motifs", plotScale="no", sortrows=T,
            clusters=1L,autoscale=F,zmin=minVal, zmax=maxVal,
            indi=F, sort_mids=T,sort_by=c(T,F,F))
if(plotPDFs==F){
  dev.off()
}


# rex Motifs
winSize=10000
rexMotifs<-import.bed(paste0(outPath,"/publicData/rexMotifs_Albritton2017_ce11.bed"))
rexMotifs<-rexMotifs[seqnames(rexMotifs)=="chrX"]
rexMotifs<-resize(rexMotifs,width=winSize,fix="center")

flankSize<-60000
smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
                  useContrasts,"_lfc.bw")
names(smcRNAseq)<-useContrasts
for(grp in useContrasts){
  rnaSeq<-import.bw(smcRNAseq[[grp]])
  cov<-coverage(rnaSeq,weight="score")
  rexMotifs<-binnedAverage(rexMotifs,cov,grp)
}

df<-data.frame(rexMotifs)
df<-pivot_longer(df,useContrasts,names_to="SMC")
p<-ggplot(df,aes(x=as.factor(score),y=value,col=SMC)) + geom_jitter(width=0.2) +
  ylab("Average RNAseq score (10kb window)") + xlab("Rex motif score") +
  ggtitle(paste0("Average RNAseq in ",winSize/1000,"kb window vs rex motif score"))

ggsave(paste0(outPath, "/plots/",outputNamePrefix,"RNAseqVmotifScore.pdf"), p,
       device="pdf",width=15,height=12,units="cm")




winSize=10000
numWins=10
bwFiles<-paste0(outPath,"/tracks/",fileNamePrefix,
                    useContrasts,"_lfc.bw")
names(bwFiles)<-useContrasts

# rex motifs
rexMotifs<-import.bed(paste0(outPath,"/publicData/rexMotifs_Albritton2017_ce11.bed"))
rexMotifs<-rexMotifs[seqnames(rexMotifs)=="chrX"]
p<-avrSignalBins(rexMotifs, bwFiles, winSize=10000,numWins=10)
p<-p+ggtitle("Average RNAseq LFC around rex motifs")
print(p)
ggsave(paste0(outPath, "/plots/",outputNamePrefix,"binnedRNAseq_nearRexMotifs.pdf"), p,
       device="pdf",width=19,height=23,units="cm")



#### rex sites
rexSites<-import.bed(paste0(outPath,"/publicData/rexSites_Albritton2017_ce11.bed"))
seqlevels(rexSites)<-seqlevels(Celegans)[1:6]
p<-avrSignalBins(rexSites, bwFiles, winSize=10000,numWins=10)
p<-p+ggtitle("Average RNAseq LFC around rex sites")
print(p)
ggsave(paste0(outPath, "/plots/",outputNamePrefix,"binnedRNAseq_nearRexSites.pdf"), p,
       device="pdf",width=19,height=23,units="cm")

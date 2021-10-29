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



pdf(file=paste0(outPath, "/plots/",outputNamePrefix,"distance2rex_",
                "_padj",padjVal, "_lfc", lfcVal,".pdf"),
    width=11,height=8,paper="a4r")

# mex Motifs
mexMotifs<-import.bed(paste0(outPath,"/publicData/rexMotifs_Albritton2017_ce11.bed"))

for(grp in useContrasts){
  dist2rex<-distanceToNearest(sigGR[[grp]],mexMotifs,ignore.strand=T)
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
  theme_bw() +xlab("Distance to mex Motif (kb)") +
  ggtitle("Distance of significantly changed genes to nearest mex motif")
print(p1)

p2<-ggplot(dfall[dfall$XvA=="X",],aes(x=dist2rex/1000,fill=upVdown))+geom_histogram() +facet_grid(upVdown~SMC) +
  theme_bw() +xlab("Distance to mex motif (kb)") +
  ggtitle("Distance of chrX significantly changed genes to nearest mex motif")
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
## seqplots for rex/mex sites------
###############-

# anchors - seqplots heatmaps ---------------------------------------------
####
## anchors
####

# mex Motifs
winSize=10000
mexMotifs<-import.bed(paste0(outPath,"/publicData/rexMotifs_Albritton2017_ce11.bed"))
mexMotifs<-mexMotifs[seqnames(mexMotifs)=="chrX"]
mexMotifs<-resize(mexMotifs,width=winSize,fix="center")
flankSize<-winSize*6

smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
                  useContrasts,"_lfc.bw")
names(smcRNAseq)<-useContrasts
if(plotPDFs==T){
  pdf(file=paste0(outPath,"/plots/",outputNamePrefix,"mexMotifs-chrX_flank",
                  flankSize/1000,"kb.pdf"), width=19,
      height=16, paper="a4")
}

if(plotPDFs==F){
  png(filename=paste0(outPath,"/plots/",outputNamePrefix,"mexMotifs-chrX_flank",
                      flankSize/1000,"kb.png"),width=19,
      height=16, units="cm", res=150)
}

############################# seqplots ---------
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


# mex Motifs
winSize=10000
mexMotifs<-import.bed(paste0(outPath,"/publicData/rexMotifs_Albritton2017_ce11.bed"))
mexMotifs<-mexMotifs[seqnames(mexMotifs)=="chrX"]
mexMotifs<-resize(mexMotifs,width=winSize,fix="center")
flankSize<-winSize*6

smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
                  useContrasts,"_lfc.bw")
names(smcRNAseq)<-useContrasts
if(plotPDFs==T){
  pdf(file=paste0(outPath,"/plots/",outputNamePrefix,"seqplots_mexMotifs-chrX_flank",
                  flankSize/1000,"kb.pdf"), width=19,
      height=16, paper="a4")
}

if(plotPDFs==F){
  png(filename=paste0(outPath,"/plots/",outputNamePrefix,"heatmap_mexMotifs-chrX_flank",
                      flankSize/1000,"kb.png"),width=19,
      height=16, units="cm", res=150)
}



#plotting heatmaps
p<-getPlotSetArray(tracks=c(smcRNAseq),
                   features=c(mexMotifs),
                   refgenome="ce11", bin=100L, xmin=flankSize,
                   xmax=flankSize, type="af",
                   xanchored=10000)


dd<-plotHeatmap(p,plotz=F)
heatmapQuantiles<-sapply(dd$HLST,quantile,c(0.05,0.95),na.rm=T)
minVal<-min(heatmapQuantiles[1,])
maxVal<-max(heatmapQuantiles[2,])
plotHeatmap(p,main="chrX mex motifs", plotScale="no", sortrows=T,
            clusters=1L,autoscale=F,zmin=minVal, zmax=maxVal,
            indi=F, sort_mids=T,sort_by=c(T,F,F),
            clspace=c("#00008B", "#FFFFE0","#8B0000"))


#plotting averages
winSize=10000
mexMotifs<-import.bed(paste0(outPath,"/publicData/rexMotifs_Albritton2017_ce11.bed"))
mexMotifs<-mexMotifs[seqnames(mexMotifs)=="chrX"]
mexMotifs<-resize(mexMotifs,width=winSize,fix="center")
flankSize<-winSize*6
smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
                  useContrasts,"_lfc.bw")
names(smcRNAseq)<-useContrasts
if(plotPDFs==F){
  dev.off()
  png(filename=paste0(outPath,"/plots/",outputNamePrefix,"lineplot_mexMotifs-chrX_flank",
                      flankSize/1000,"kb.png"),width=19,
      height=16, units="cm", res=150)
}

p<-getPlotSetArray(tracks=c(smcRNAseq),
                   features=c(mexMotifs),
                   refgenome="ce11", bin=winSize/10, xmin=flankSize,
                   xmax=flankSize, type="af",
                   xanchored=winSize)
dd<-plotHeatmap(p,plotz=F)
heatmapQuantiles<-sapply(dd$HLST,quantile,c(0.05,0.95),na.rm=T)
minVal<-min(heatmapQuantiles[1,])
maxVal<-max(heatmapQuantiles[2,])
plotAverage(p,main="chrX mex motifs",plotScale="",ylim=c(minVal,maxVal),
            error.estimates=ifelse(length(useContrasts>3),F,T))


dev.off()


# rex sites  -----
winSize=1000
rexSites<-import.bed(paste0(outPath,"/publicData/rexSites_Albritton2017_ce11.bed"))
rexSites<-rexSites[seqnames(rexSites)=="chrX"]
rexSites<-resize(rexSites,width=winSize,fix="center")
flankSize<-winSize*6

smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
                  useContrasts,"_lfc.bw")
names(smcRNAseq)<-useContrasts
if(plotPDFs==T){
  pdf(file=paste0(outPath,"/plots/",outputNamePrefix,"seqplots_rexSites-chrX_flank",
                  flankSize/1000,"kb.pdf"), width=19,
      height=16, paper="a4")
}

if(plotPDFs==F){
  png(filename=paste0(outPath,"/plots/",outputNamePrefix,"heatmap_rexSites-chrX_flank",
                      flankSize/1000,"kb.png"),width=19,
      height=16, units="cm", res=150)
}



#plotting heatmaps
p<-getPlotSetArray(tracks=c(smcRNAseq),
                   features=c(rexSites),
                   refgenome="ce11", bin=100L, xmin=flankSize,
                   xmax=flankSize, type="af",
                   xanchored=10000)


dd<-plotHeatmap(p,plotz=F)
heatmapQuantiles<-sapply(dd$HLST,quantile,c(0.05,0.95),na.rm=T)
minVal<-min(heatmapQuantiles[1,])
maxVal<-max(heatmapQuantiles[2,])
plotHeatmap(p,main="chrX rex sites", plotScale="no", sortrows=T,
            clusters=1L,autoscale=F,zmin=minVal, zmax=maxVal,
            indi=F, sort_mids=T,sort_by=c(T,F,F),
            clspace=c("#00008B", "#FFFFE0","#8B0000"))

#plotting averages
winSize=10000
rexSites<-import.bed(paste0(outPath,"/publicData/rexSites_Albritton2017_ce11.bed"))
rexSites<-rexSites[seqnames(rexSites)=="chrX"]
rexSites<-resize(rexSites,width=winSize,fix="center")
flankSize<-winSize*6
smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
                  useContrasts,"_lfc.bw")
names(smcRNAseq)<-useContrasts
if(plotPDFs==F){
  dev.off()
  png(filename=paste0(outPath,"/plots/",outputNamePrefix,"lineplot_rexSites-chrX_flank",
                      flankSize/1000,"kb.png"),width=19,
      height=16, units="cm", res=150)
}

p<-getPlotSetArray(tracks=c(smcRNAseq),
                   features=c(rexSites),
                   refgenome="ce11", bin=winSize/10, xmin=flankSize,
                   xmax=flankSize, type="af",
                   xanchored=winSize)
p$data$feature_1
dd<-plotHeatmap(p,plotz=F)
heatmapQuantiles<-sapply(dd$HLST,quantile,c(0.05,0.95),na.rm=T)
minVal<-min(heatmapQuantiles[1,])
maxVal<-max(heatmapQuantiles[2,])
plotAverage(p,main="chrX rex sites",plotScale="",ylim=c(minVal,maxVal),
            error.estimates=ifelse(length(useContrasts>3),F,T))


dev.off()





#############################-
## Motif score vs average RNAseq score-----
#############################-
# mex Motifs
winSize=1000
mexMotifs<-import.bed(paste0(outPath,"/publicData/rexMotifs_Albritton2017_ce11.bed"))
mexMotifs<-mexMotifs[seqnames(mexMotifs)=="chrX"]
mexMotifs<-resize(mexMotifs,width=winSize,fix="center")

#flankSize<-60000
smcRNAseq<-paste0(outPath,"/tracks/",fileNamePrefix,
                  useContrasts,"_lfc.bw")
names(smcRNAseq)<-useContrasts
for(grp in useContrasts){
  rnaSeq<-import.bw(smcRNAseq[[grp]])
  cov<-coverage(rnaSeq,weight="score")
  mexMotifs<-binnedAverage(mexMotifs,cov,grp)
}

df<-data.frame(mexMotifs)
colnames(df)<-c(colnames(df)[1:5],colnames(mcols(mexMotifs)))
df<-pivot_longer(df,useContrasts,names_to="SMC")
p<-ggplot(df,aes(x=as.factor(score),y=value,col=SMC)) + geom_jitter(width=0.2) +
  ylab(paste0("Average RNAseq score (",winSize/1000,"kb window)")) + xlab("Mex motif score") +
  ggtitle(paste0("Average RNAseq in ",winSize/1000,"kb window vs mex motif score"))

ggsave(paste0(outPath, "/plots/",outputNamePrefix,"RNAseqVmotifScore.pdf"), p,
       device="pdf",width=15,height=12,units="cm")




winSize=10000
numWins=10
bwFiles<-paste0(outPath,"/tracks/",fileNamePrefix,
                    useContrasts,"_lfc.bw")
names(bwFiles)<-useContrasts

# mex motifs
mexMotifs<-import.bed(paste0(outPath,"/publicData/rexMotifs_Albritton2017_ce11.bed"))
mexMotifs<-mexMotifs[seqnames(mexMotifs)=="chrX"]
p<-avrSignalBins(mexMotifs, bwFiles, winSize=winSize,numWins=numWins)
p<-p+ggtitle("Average RNAseq LFC around mex motifs")
print(p)
ggsave(paste0(outPath, "/plots/",outputNamePrefix,"binnedRNAseq_nearmexMotifs.pdf"), p,
       device="pdf",width=19,height=23,units="cm")



#### rex sites
rexSites<-import.bed(paste0(outPath,"/publicData/rexSites_Albritton2017_ce11.bed"))
seqlevels(rexSites)<-seqlevels(Celegans)[1:6]
p<-avrSignalBins(rexSites, bwFiles, winSize=winSize,numWins=numWins)
p<-p+ggtitle("Average RNAseq LFC around rex sites")
print(p)
ggsave(paste0(outPath, "/plots/",outputNamePrefix,"binnedRNAseq_nearRexSites.pdf"), p,
       device="pdf",width=19,height=23,units="cm")


#################3
#################
states<-import.bed("./publicData/chromStates_L3_Evans2016_ce11.bed")

stateClrs<-c("#fe0003","#f59745","#008100","#74943a",
             "#c4d69c","#05ff7f","#ceff65","#fd0082",
             "#ff70d0","#ffb6c4","#7f00ff","#1900fe",
             "#528dd4","#814006","#b8cde4","#808080",
             "#dcdac3","#c4bc97","#938b54","#141414")
names(stateClrs)<-c("Promoter","5pProx&gene","TxnElongationI","TxnElongationII",
                    "TxnElongationIII","EnhancerI","EnhancerII","EnhancerIII",
                    "EnhancerIV","EnahncerV","LowExpn","H3K9me2","Repeats",
                    "TissueSpecific_high","Tissuespecific_low","TissueSpecific_moderate",
                    "Repeats&intergenic","H3K27me3_I","H3K27me3_II","H3K9me3")

if(!file.exists("./publicData/chromStates_L3_Evans2016_ce11_rgb.bed")){
  stateRGB<-apply(col2rgb(stateClrs),2,paste,collapse=",")
  states2bed<-states
  md<-data.frame(name=states2bed$score, score=states2bed$score,
                 strand=".",
                 thickStart=start(states), thickEnd=end(states),
                 itemRGB=stateRGB[states$score],blockCount="1",
                 blockSizes=width(states),blockStarts=start(states))
  mcols(states2bed)<-md
  trackLine<-'track name="chromStates" description="Evans et al. 2016" visibility=1 itemRgb="On"\n'
  export.bed(states2bed,"./publicData/chromStates_L3_Evans2016_ce11_rgb.bed")
  states2bed<-read.delim("./publicData/chromStates_L3_Evans2016_ce11_rgb.bed",header=F)
  states2bed<-cbind(states2bed,md[,c("itemRGB")])
  write.table(states2bed,
              file="./publicData/chromStates_L3_Evans2016_ce11_rgb.bed",
              sep="\t",row.names=F,col.names=F,quote=F)
  # add trackline (but not necessary)
  fConn <- file("./publicData/chromStates_L3_Evans2016_ce11_rgb.bed", 'r+')
  Lines <- readLines(fConn)
  writeLines(c(trackLine, Lines), con = fConn)
  close(fConn)
}




######### mex motif score vs distance to nearest gene
md<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))
mexMotifs<-import.bed(paste0(outPath,"/publicData/rexMotifs_Albritton2017_ce11.bed"))
rexSites<-import.bed(paste0(outPath,"/publicData/rexSites_Albritton2017_ce11.bed"))
seqlevels(mexMotifs)<-seqlevels(Celegans)
seqlevels(rexSites)<-seqlevels(Celegans)
mexMotifs<-mexMotifs[seqnames(mexMotifs)=="chrX"]

# find mex in rex
ol<-findOverlaps(mexMotifs,rexSites)
mexMotifs$mexInRex<-F
mexMotifs$mexInRex[unique(queryHits(ol))]<-T


# find distance to nearest
dtn<-distanceToNearest(mexMotifs,md,ignore.strand=T)
mexMotifs$dist_anyStrand<-mcols(dtn)$distance
dtn<-distanceToNearest(mexMotifs,md,ignore.strand=F)
mexMotifs$dist_sameStrand<-mcols(dtn)$distance

# find expression of nearest
sigGR<-list()
for (grp in useContrasts){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],
                         "_DESeq2_fullResults_p",padjVal,".rds"))
  sigGR[[grp]]<-GRanges(salmon)
}

#anyStrand
mexMotifs$nearestGene_anyStrand<-md$wormbaseID[nearest(mexMotifs,md,ignore.strand=T)]
mexMotifs$lfcNearest_anyStrand<-0
idx1<-mexMotifs$nearestGene_anyStrand %in% sigGR[["dpy26"]]$wormbaseID
idx2<-match(mexMotifs$nearestGene_anyStrand[idx1],sigGR[["dpy26"]]$wormbaseID)
mexMotifs$lfcNearest_anyStrand[idx1]<-sigGR[["dpy26"]]$log2FoldChange[idx2]

#sameStrand
mexMotifs$nearestGene_sameStrand<-md$wormbaseID[nearest(mexMotifs,md,ignore.strand=F)]
mexMotifs$lfcNearest_sameStrand<-0
idx1<-mexMotifs$nearestGene_sameStrand %in% sigGR[["dpy26"]]$wormbaseID
idx2<-match(mexMotifs$nearestGene_sameStrand[idx1],sigGR[["dpy26"]]$wormbaseID)
mexMotifs$lfcNearest_sameStrand[idx1]<-sigGR[["dpy26"]]$log2FoldChange[idx2]


# find chromatin states
ol<-findOverlaps(mexMotifs,states)
mexMotifs$chromState<-NA
mexMotifs$chromState[queryHits(ol)]<-states$score[subjectHits(ol)]

df<-data.frame(mexMotifs)
df$chromState<-factor(df$chromState)
levels(df$chromState)<-c(names(stateClrs))

pdf(paste0(outPath,"/plots/",outputNamePrefix,"MexMotifScoreVlfc.pdf"),width=11,height=8, paper="a4r")

p1<-ggplot(df,aes(x=score,y=lfcNearest_anyStrand,shape=mexInRex,col=chromState,
                  size=dist_anyStrand)) + geom_point() +
  scale_color_manual(values=stateClrs,na.value="white")+
  ggtitle("Mex motif score vs LFC of nearest gene (both strands)")
print(p1)
p2<-ggplot(df,aes(x=score,y=lfcNearest_sameStrand,shape=mexInRex,col=chromState,
                  size=dist_sameStrand)) + geom_point()+
  scale_color_manual(values=c(stateClrs),na.value="white")+
  ggtitle("Mex motif score vs LFC of nearest gene (same strand)")
print(p2)
dev.off()

pdf(paste0(outPath,"/plots/",outputNamePrefix,"MexMotifScoreVdistance.pdf"),width=11,height=8, paper="a4r")

p1<-ggplot(df,aes(x=score,y=dist_anyStrand,shape=mexInRex,col=chromState,
                  size=lfcNearest_anyStrand)) + geom_point() +
  scale_color_manual(values=stateClrs,na.value="white")+
  ggtitle("Mex motif score vs distance to nearest gene (both strands)")
print(p1)
p2<-ggplot(df,aes(x=score,y=dist_sameStrand,shape=mexInRex,col=chromState,
                  size=lfcNearest_sameStrand)) + geom_point()+
  scale_color_manual(values=c(stateClrs),na.value="white")+
  ggtitle("Mex motif score vs distance to nearest gene (same strand)")
print(p2)
dev.off()



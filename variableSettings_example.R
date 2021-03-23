library(GenomicRanges)

plotPDFs=T
padjVal=0.05
lfcVal=0.5
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"/salmon_")
filterPrefix=paste0("p",padjVal,"_lfc",lfcVal,"/prefiltCyc2xChrAX_")

outPath="."
genomeVer="WS275"
genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)

remakeFiles=F # remake publicData files?
combineChrAX=F # artificially combine chrA and X from different datasets
filterData=T
if(filterData){
  oscillating<-read.delim(paste0(outPath,"/publicData/oscillatingGenes.tsv"), header=T,
                          stringsAsFactors=F) #3739
  latorre<-read.delim(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv")) #3235
  #hsUP<-readRDS(file=paste0(outPath,"/publicData/hsUp_garrigues2019.rds")) #1680
  #hsDOWN<-readRDS(file=paste0(outPath,"/publicData/hsDown_garrigues2019.rds")) #455
  #toFilter<-unique(c(oscillating$wormbaseID, latorre$wormbaseID,
  #hsUP$wormbaseID, hsDOWN$wormbaseID))
  #md<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))
  #chrXidx<-as.vector(seqnames(md)=="chrX")
  #length(md$wormbaseID[chrXidx]) #5821
  toFilter<-unique(c(oscillating$wormbaseID,latorre$wormbaseID))
  print(paste0("filtering ", length(toFilter), " genes"))
  #4522 genes osc+latorre
  #9541 genes osc+latorre+chrX
  #6101 genes osc+latorre+hs
}



plotPDFs=F
padjVal=0.05
lfcVal=0.5
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"/salmon_")
filterPrefix=paste0("p",padjVal,"_lfc",lfcVal,"/noOsc_")

filterData=F
if(filterData){
  oscillating<-read.delim(paste0(outPath,"/oscillatingGenes.tsv"), header=T,
                          stringsAsFactors=F) #3739
  latorre<-read.delim(paste0(outPath,"/oscillatingGenes_latorre.tsv")) #3235
  #hsUP<-readRDS(file="hsUp_garrigues2019.rds") #1680
  #hsDOWN<-readRDS(file="hsDown_garrigues2019.rds") #455
  #toFilter<-unique(c(oscillating$WB_ID, latorre$wormbaseID,
  #hsUP$WormBase.ID, hsDOWN$WormBase.ID))
  toFilter<-unique(c(oscillating$WB_ID,latorre$wormbaseID))
  #4522 genes osc+latorre
  #6101 genes osc+latorre+hs
}

outPath="."
genomeVer="WS275"
genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)

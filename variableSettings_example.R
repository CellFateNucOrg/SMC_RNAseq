library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)

plotPDFs=F
padjVal=0.05
lfcVal=0.5
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"/preFiltChrAX_")
filterPrefix=paste0("p",padjVal,"_lfc",lfcVal,"/preFiltChrAX_")

outPath="."
genomeVer="WS275"
genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)

wbseqinfo<-seqinfo(Celegans)
seqnames(wbseqinfo)<-c(gsub("chr","",seqnames(Celegans)))
seqnames(wbseqinfo)<-c(gsub("^M$","MtDNA",seqnames(wbseqinfo)))
genome(wbseqinfo)<-genomeVer
ce11seqinfo<-seqinfo(Celegans)

remakeFiles=F # remake publicData files?
combineChrAX=T # artificially combine chrA and X from different datasets?
if(combineChrAX){
  chrXprefix<-paste0("/../SMC_RNAseq_prefilt/rds/p",padjVal,"_lfc",lfcVal,"/preFilt_")
  chrAprefix<-paste0("/../SMC_RNAseq_prefiltChrA/rds/p",padjVal,"_lfc",lfcVal,"/preFiltChrA_")
}

filterData=F
if(filterData){
    #oscillating<-read.delim(paste0(outPath,"/publicData/oscillatingGenes.tsv"), header=T,
    #                        stringsAsFactors=F) #3739
    #latorre<-read.delim(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv")) #3235
    #hsUP<-readRDS(file=paste0(outPath,"/publicData/hsUp_garrigues2019.rds")) #1680
    #hsDOWN<-readRDS(file=paste0(outPath,"/publicData/hsDown_garrigues2019.rds")) #455
    #toFilter<-unique(c(oscillating$wormbaseID, latorre$wormbaseID,
    #hsUP$wormbaseID, hsDOWN$wormbaseID))
  #md<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))
  #chrXidx<-as.vector(seqnames(md)=="chrX")
  #length(md$wormbaseID[chrXidx]) #5821
  #toFilter<-unique(c(md$wormbaseID[chrXidx]))
  #print(paste0("filtering ", length(toFilter), " genes"))
  #4522 genes osc+latorre
  #9541 genes osc+latorre+chrX
  #6101 genes osc+latorre+hs
}


strainLevels<-c("366","382","775","784")
varOIlevels<-c("wt","dpy26cs","kle2cs","scc1cs")
#strainLevels<-c("366","821","823")
#varOIlevels<-c("wt","dpy26cs_sdc3deg","TIR")
varOI<-"SMC"

fileList<-read.table(paste0(outPath,"/fastqList.txt"),stringsAsFactors=F,header=T)

sampleNames<-paste(fileList$sampleName, fileList$repeatNum, fileList$laneNum, sep="_")

fileNames<-paste0(outPath,"/salmon/mRNA/",sampleNames,"/quant.sf")
sampleTable<-data.frame(fileName=fileNames,sampleName=sampleNames,stringsAsFactors=F)

# extract the technical replicate variable
sampleTable$replicate=factor(fileList$repeatNum)
sampleTable$lane=factor(fileList$laneNum)

# extract the strain variable
sampleTable$strain<-factor(as.character(fileList$sampleName),levels=strainLevels)
sampleTable[,varOI]<-sampleTable$strain
levels(sampleTable[,varOI])<-varOIlevels

controlGrp<-levels(sampleTable[,varOI])[1] # control group
groupsOI<-levels(sampleTable[,varOI])[-1] # groups of interest to contrast to control



library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)

plotPDFs=F
padjVal=0.05
lfcVal=0.5
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"/salmon_")
filterPrefix=paste0("p",padjVal,"_lfc",lfcVal,"/filtCycChrAX_")

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
  chrXprefix<-paste0("/rds/p",padjVal,"_lfc",lfcVal,"_filtCyc/filtCyc_")
  chrAprefix<-paste0("/rds/p",padjVal,"_lfc",lfcVal,"_filtCycChrA/filtCycChrA_")
}

filterData=T
if(filterData){
  oscillating<-read.delim(paste0(outPath,"/publicData/oscillatingGenes.tsv"), header=T,
                            stringsAsFactors=F) #3739
  latorre<-read.delim(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv")) #3235
  toFilter<-unique(c(oscillating$wormbaseID, latorre$wormbaseID))
  #hsUP<-readRDS(file=paste0(outPath,"/publicData/hsUp_garrigues2019.rds")) #1680
  #hsDOWN<-readRDS(file=paste0(outPath,"/publicData/hsDown_garrigues2019.rds")) #455
  #toFilter<-unique(c(oscillating$wormbaseID, latorre$wormbaseID,
  #hsUP$wormbaseID, hsDOWN$wormbaseID))
  #md<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))
  #chrXidx<-as.vector(seqnames(md)=="chrX")
  #length(md$wormbaseID[chrXidx]) #5821
  #toFilter<-unique(c(oscillating$wormbaseID, latorre$wormbaseID, md$wormbaseID[chrXidx]))
  #print(paste0("filtering ", length(toFilter), " genes"))
  #4522 genes osc+latorre
  #9541 genes osc+latorre+chrX
  #6101 genes osc+latorre+hs
}


#strainLevels<-c("366","382","775","784")
#varOIlevels<-c("wt","dpy26cs","kle2cs","scc1cs")
#varOIlevels<-c("wt","dpy26cs_sdc3deg","TIR")
strainLevels<-c("366","382","821","822","823")
varOIlevels<-c("wt.wt.wt.0mM","wt.wt.wt.1mM","wt.TIR1.wt.1mM" ,"wt.TIR1.sdc3deg.0mM","wt.TIR1.sdc3deg.1mM","dpy26cs.wt.wt.0mM","dpy26cs.TIR1.sdc3deg.1mM")
varOI<-"SMC"
#strainLevels<-c("366","822","823")
#varOIlevels<-c("wt.0mM.wt","wt.1mM.wt","TIR1.1mM.wt",
#               "TIR1.0mM.sdc3deg","TIR1.1mM.sdc3deg")
#varOI<-"sdc3deg"



fileList<-read.table(paste0(outPath,"/fastqList.txt"),stringsAsFactors=F,header=T)

#fileList<-fileList[fileList$sampleName %in% strainLevels,]
#fileList<-fileList[fileList$date != 20201023,]

sampleNames<-paste(fileList$sampleName, fileList$repeatNum, fileList$laneNum, sep="_")

fileNames<-paste0(outPath,"/salmon/mRNA/",sampleNames,"/quant.sf")
sampleTable<-data.frame(fileName=fileNames,sampleName=sampleNames,stringsAsFactors=F)
sampleTable$strain<-factor(as.character(fileList$sampleName),levels=strainLevels)

# extract the technical replicate variable
sampleTable$date<-factor(fileList$date)
sampleTable$batch<-factor(fileList$repeatNum)
sampleTable$replicate<-gsub("^.?","",fileList$repeatNum)
sampleTable$lane<-factor(gsub("a","",fileList$laneNum))
sampleTable$dpy26<-"wt"
sampleTable$dpy26[sampleTable$strain %in% c("821","382")]<-"dpy26cs"
sampleTable$dpy26<-factor(sampleTable$dpy26, levels=c("wt","dpy26cs"))

sampleTable$TIR1<-"TIR1"
sampleTable$TIR1[sampleTable$strain %in% c("366","382")]<-"wt"
sampleTable$TIR1<-factor(sampleTable$TIR1, levels=c("wt","TIR1"))

sampleTable$sdc3<-"wt"
sampleTable$sdc3[sampleTable$strain %in% c("821","822")]<-"sdc3deg"
sampleTable$sdc3<-factor(sampleTable$sdc3, levels=c("wt","sdc3deg"))

sampleTable$auxin<-"1mM"
sampleTable$auxin[grep("C|B",sampleTable$batch)]<-"0mM"
sampleTable$auxin<-factor(sampleTable$auxin,levels=c("0mM","1mM"))

sampleTable$sdc3deg<-factor(paste(sampleTable$TIR1, sampleTable$auxin ,sampleTable$sdc3,sep="."),
                            levels=c("wt.0mM.wt","wt.1mM.wt","TIR1.1mM.wt",
                                     "TIR1.0mM.sdc3deg","TIR1.1mM.sdc3deg"))

sampleTable$auxTIR<-factor(paste(sampleTable$auxin,sampleTable$TIR1,sep="."),
                           levels=c("0mM.wt","1mM.wt","0mM.TIR1","1mM.TIR1"))


# extract the strain variable
sampleTable[,varOI]<-factor(paste(sampleTable$dpy26, sampleTable$TIR1, sampleTable$sdc3, sampleTable$auxin, sep = "."),levels=varOIlevels)

controlGrp<-levels(sampleTable[,varOI])[1] # control group
groupsOI<-levels(sampleTable[,varOI])[-1] # groups of interest to contrast to control

###################
## building model and contrasts of interest
###################

# proposed model
modelTxt<-"~replicate+lane+date+SMC"
# to use as a model, must be converted to formula, i.e. formula(modelTxt)

# to get names of contrasts of interest that are already in the model:
contrastsOI<-as.list(paste0(varOI,"_",groupsOI,"_vs_",controlGrp))
names(contrastsOI)<-contrastsOI
contrastNames<-as.list(names(contrastsOI))
names(contrastNames)<-groupsOI
useContrasts<-groupsOI

##### more complex contrasts: ##############
# for advanced contrast matrix usage. see:
# https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
# Create the model.matrix (this is similar but not quite the same as betting it from dds)
advancedContrasts=T
if(advancedContrasts){
  contrastsOI<-list()
  modMat<-model.matrix(formula(modelTxt),data=sampleTable)
  colnames(modMat)
  # set control variables columns from matrix to 0
  modMat[,grep(paste0("^",varOI),colnames(modMat),invert=T)]<-0

  # get coefficients for existing subsets in data table
  wt_wt_wt_0mM<-colMeans(modMat[sampleTable$SMC == "wt.wt.wt.0mM", ])
  wt_wt_wt_1mM<-colMeans(modMat[sampleTable$SMC == "wt.wt.wt.1mM", ])
  wt_TIR1_wt_1mM<-colMeans(modMat[sampleTable$SMC == "wt.TIR1.wt.1mM", ])
  wt_TIR1_sdc3deg_0mM<-colMeans(modMat[sampleTable$SMC == "wt.TIR1.sdc3deg.0mM", ])
  wt_TIR1_sdc3deg_1mM<-colMeans(modMat[sampleTable$SMC == "wt.TIR1.sdc3deg.1mM", ])
  dpy26cs_wt_wt_0mM<-colMeans(modMat[sampleTable$SMC == "dpy26cs.wt.wt.0mM", ])
  dpy26cs_TIR1_sdc3deg_1mM<-colMeans(modMat[sampleTable$SMC == "dpy26cs.TIR1.sdc3deg.1mM", ])

  # add special contrasts by subtracting relevant coefficient vectors
  #contratsOI<-list()
  contrastsOI[["wt.wt.wt.X_1mM_vs_0mM"]]<-wt_wt_wt_1mM #aux
  contrastsOI[["wt.TIR1.sdc3deg.X_1mM_vs_0mM"]]<-wt_TIR1_sdc3deg_1mM - wt_TIR1_sdc3deg_0mM #aux_sdc3BG
  contrastsOI[["wt.X.wt.1mM_TIR1_vs_wt"]]<-wt_TIR1_wt_1mM - wt_wt_wt_1mM #TIR1
  contrastsOI[["wt.X.wt.X_TIR11mM_vs_wt0mM"]]<-wt_TIR1_wt_1mM - wt_wt_wt_0mM #TIRaux
  contrastsOI[["wt.TIR1.X.1mM_sdc3deg_vs_wt"]]<-wt_TIR1_sdc3deg_1mM - wt_TIR1_wt_1mM #sdc3
  contrastsOI[["X.wt.wt.0mM_dpy26cs_vs_wt"]]<-dpy26cs_wt_wt_0mM - wt_wt_wt_0mM #dpy26
  contrastsOI[["X.TIR1.sdc3deg.1mM_dpy26cs_vs_wt"]]<-dpy26cs_TIR1_sdc3deg_1mM - wt_TIR1_sdc3deg_1mM #dpy26_sdc3BG
  contrastsOI[["dpy26cs.X.X.X_TIR1sdc3deg1mM_vs_wtwt0mM"]]<-dpy26cs_TIR1_sdc3deg_1mM - wt_TIR1_wt_1mM - dpy26cs_wt_wt_0mM #sdc3_dpy26BG
  contrastsOI[["X.TIR1.X.1mM_dpy26cssdc3deg_vs_wtwt"]]<-dpy26cs_TIR1_sdc3deg_1mM - wt_TIR1_wt_1mM #dpy26sdc3

  contrastNames<-as.list(names(contrastsOI))
  names(contrastNames)<-c("aux","aux_sdc3BG","TIR1","TIR1aux","sdc3","dpy26","dpy26_sdc3BG","sdc3_dpy26BG","dpy26sdc3")

  useContrasts<-c("aux_sdc3BG","sdc3","dpy26")
}

advancedContrasts=F
if(advancedContrasts){
  modMat<-model.matrix(formula(modelTxt),data=sampleTable)
  colnames(modMat)
  # set control variables columns from matrix to 0
  modMat[,grep(paste0("^",varOI),colnames(modMat),invert=T)]<-0

  # get coefficients for existing subsets in data table
  #wt.0mM.wt<-colMeans(modMat[sampleTable[,varOI] == "wt.0mM.wt", ])
  wt.1mM.wt<-colMeans(modMat[sampleTable[,varOI] == "wt.1mM.wt", ])
  TIR1.1mM.wt<-colMeans(modMat[sampleTable[,varOI] == "TIR1.1mM.wt", ])
  TIR1.0mM.sdc3deg<-colMeans(modMat[sampleTable[,varOI] == "TIR1.0mM.sdc3deg", ])
  TIR1.1mM.sdc3deg<-colMeans(modMat[sampleTable[,varOI] == "TIR1.1mM.sdc3deg", ])

  # add special contrasts by subtracting relevant coefficient vectors
  #contratsOI<-list()
  contrastsOI[["wt.X.wt_1mM_vs_0mM"]]<-wt.1mM.wt #aux
  contrastsOI[["TIR1.X.sdc3deg_1mM_vs_0mM"]]<-TIR1.1mM.sdc3deg - TIR1.0mM.sdc3deg #aux_sdc3BG
  contrastsOI[["X.1mM.wt_TIR1_vs_wt"]]<-TIR1.1mM.wt - wt.1mM.wt #TIR1
  contrastsOI[["X.X.wt_TIR11mM_vs_wt0mM"]]<-TIR1.1mM.wt #TIR1aux
  contrastsOI[["TIR1.1mM.X_sdc3deg_vs_wt"]]<-TIR1.1mM.sdc3deg - TIR1.1mM.wt #sdc3

  contrastNames<-as.list(names(contrastsOI))
  names(contrastNames)<-c("aux","aux_sdc3BG","TIR1","TIR1aux","sdc3")

  useContrasts<-c("aux_sdc3BG","sdc3")
}

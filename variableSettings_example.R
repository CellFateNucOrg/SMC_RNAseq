library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library(dplyr)

plotPDFs=F
padjVal=0.05
lfcVal=0.5
fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"/no775B3_")
outPath="."
genomeVer="WS275"
genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)
fileList<-read.table(paste0(outPath,"/fastqList_no775B3.txt"),stringsAsFactors=F,header=T)

rnaType="mRNA"
remakeFiles=F # remake publicData files?
combineChrAX=F # artificially combine chrA and X from different datasets?
filterData=T # filter by certain gene lists
#filterBy=c("chrX") #names in filterList of gene lists to use
filterBy=c("Cycling_Meeuse","Cycling_Latorre")
#filterBy=c("Cycling_Meeuse","Cycling_Latorre","chrX")
customNameTxt=ifelse(rnaType=="mRNA","",rnaType) #some other text you want to add to the filename prefix


# get ce11 genome
wbseqinfo<-seqinfo(Celegans)
seqnames(wbseqinfo)<-c(gsub("chr","",seqnames(Celegans)))
seqnames(wbseqinfo)<-c(gsub("^M$","MtDNA",seqnames(wbseqinfo)))
genome(wbseqinfo)<-genomeVer
ce11seqinfo<-seqinfo(Celegans)

###### don't change these variables
filterName=NULL
filterPrefix=NULL
chrXprefix=NULL
chrAprefix=NULL
filterList<-list()
if(filterData | combineChrAX){
  filterList[["Cycling_Meeuse"]]<-read.delim(paste0(outPath,"/publicData/oscillatingGenes.tsv"), header=T, stringsAsFactors=F)$wormbaseID #3739
  filterList[["Cycling_Latorre"]]<-read.delim(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv"))$wormbaseID #3235
  #hsUP<-readRDS(file=paste0(outPath,"/publicData/hsUp_garrigues2019.rds")) #1680
  #hsDOWN<-readRDS(file=paste0(outPath,"/publicData/hsDown_garrigues2019.rds")) #455
  #toFilter<-unique(c(oscillating$wormbaseID, latorre$wormbaseID,
  #hsUP$wormbaseID, hsDOWN$wormbaseID))
  md<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))
  chrXidx<-as.vector(seqnames(md)=="chrX")
  #length(md$wormbaseID[chrXidx]) #5821
  filterList[["chrX"]]<-md$wormbaseID[chrXidx]
  #print(paste0("filtering ", length(toFilter), " genes"))
  #4522 genes osc+latorre
  #9541 genes osc+latorre+chrX
  #6101 genes osc+latorre+hs
  lapply(filterList,length)
  toFilter<-unique(unlist(filterList[filterBy]))
  length(toFilter) #9603
  filterName<-paste0("filt",
                       ifelse(any(grepl("Cycling",filterBy)),"Cyc",""),
                       ifelse(combineChrAX,"ChrAX",""),
                       ifelse(any(grepl("chrX",filterBy)) & !combineChrAX,"ChrA",""),
                       ifelse(customNameTxt!="",paste0("_",customNameTxt),""))
  fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"_",filterName,"/",
                        ifelse(customNameTxt!="",customNameTxt,"salmon"),"_")
  filterPrefix=paste0("p",padjVal,"_lfc",lfcVal,"_",filterName,"/",
                      filterName,"_")
}

if(combineChrAX){
  chrXprefix<-paste0("/rds/",gsub("ChrAX","",filterPrefix))
  chrXprefix<-ifelse(grepl("_filt\\/filt_$",chrXprefix),paste0("/rds/","p",padjVal,"_lfc",lfcVal,"/no775B3_"),chrXprefix)
  chrAprefix<-paste0("/rds/",gsub("ChrAX","ChrA",filterPrefix))
  fileNamePrefix<-filterPrefix
}


########## getting variables for model-----

strainLevels<-c("366","382","821","822","823","775","784","828","844")
varOIlevels<-c("wt.wt.wt.0mM","wt.wt.wt.1mM","wt.TIR1.wt.1mM" ,"wt.TIR1.sdc3deg.0mM","wt.TIR1.sdc3deg.1mM","dpy26cs.wt.wt.0mM","dpy26cs.TIR1.sdc3deg.1mM","kle2cs.wt.wt.0mM","scc1cs.wt.wt.0mM","coh1cs.wt.wt.0mM","scc1coh1cs.wt.wt.0mM")
varOI<-"SMC"


sampleNames<-paste(fileList$sampleName, fileList$repeatNum, fileList$laneNum, sep="_")

fileNames<-paste0(outPath,"/salmon/",rnaType,"/",sampleNames,"/quant.sf")
sampleTable<-data.frame(fileName=fileNames,sampleName=sampleNames,stringsAsFactors=F)
sampleTable$strain<-factor(as.character(fileList$sampleName),levels=strainLevels)

# extract the technical replicate variable
sampleTable$date<-factor(fileList$date)
sampleTable$batch<-factor(fileList$repeatNum)
sampleTable$replicate<-gsub("^.?","",fileList$repeatNum)
sampleTable$replicateBatch<-gsub(".?$","",fileList$repeatNum)
sampleTable$lane<-factor(gsub("a","",fileList$laneNum))

#kleisin
sampleTable$kleisin<-"wt"
sampleTable$kleisin[sampleTable$strain %in% c("821","382")]<-"dpy26cs"
sampleTable$kleisin[sampleTable$strain %in% c("775")]<-"kle2cs"
sampleTable$kleisin[sampleTable$strain %in% c("784")]<-"scc1cs"
sampleTable$kleisin[sampleTable$strain %in% c("828")]<-"coh1cs"
sampleTable$kleisin[sampleTable$strain %in% c("844")]<-"scc1coh1cs"
sampleTable$kleisin<-factor(sampleTable$kleisin, levels=c("wt","dpy26cs",
                                                        "kle2cs","scc1cs",
                                                        "coh1cs","scc1coh1cs"))
#TIR1
sampleTable$TIR1<-"TIR1"
sampleTable$TIR1[sampleTable$strain %in% c("366","382","775","784","828","844")]<-"wt"
sampleTable$TIR1<-factor(sampleTable$TIR1, levels=c("wt","TIR1"))

#sdc3
sampleTable$sdc3<-"wt"
sampleTable$sdc3[sampleTable$strain %in% c("821","822")]<-"sdc3deg"
sampleTable$sdc3<-factor(sampleTable$sdc3, levels=c("wt","sdc3deg"))

#auxin
sampleTable$auxin<-"1mM"
sampleTable$auxin[grep("C|B",sampleTable$batch)]<-"0mM"
sampleTable$auxin<-factor(sampleTable$auxin,levels=c("0mM","1mM"))


# extract the strain variable
sampleTable[,varOI]<-factor(paste(sampleTable$kleisin, sampleTable$TIR1, sampleTable$sdc3, sampleTable$auxin, sep = "."),levels=varOIlevels)

controlGrp<-levels(sampleTable[,varOI])[1] # control group
groupsOI<-levels(sampleTable[,varOI])[-1] # groups of interest to contrast to control

sampleTable<-sampleTable %>% dplyr::group_by(strain,replicateBatch,lane) %>% dplyr::mutate(count=n()) %>% dplyr::mutate(newRep=1:count)
sampleTable$count<-factor(sampleTable$count)
sampleTable$newRep<-factor(sampleTable$newRep)

###################-
## building model and contrasts of interest------
###################-

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
# Create the model.matrix (this is similar but not quite the same as getting it from dds)
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
  kle2cs_wt_wt_0mM<-colMeans(modMat[sampleTable$SMC == "kle2cs.wt.wt.0mM", ])
  scc1cs_wt_wt_0mM<-colMeans(modMat[sampleTable$SMC == "scc1cs.wt.wt.0mM", ])
  coh1cs_wt_wt_0mM<-colMeans(modMat[sampleTable$SMC == "coh1cs.wt.wt.0mM", ])
  scc1coh1cs_wt_wt_0mM<-colMeans(modMat[sampleTable$SMC == "scc1coh1cs.wt.wt.0mM", ])
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
  contrastsOI[["X.TIR1.sdc3deg.X_dpy26csaux_vs_wt0mM"]]<-dpy26cs_TIR1_sdc3deg_1mM - wt_TIR1_sdc3deg_0mM #aux_dpy26sdc3
  contrastsOI[["X.wt.wt.0mM_kle2cs_vs_wt"]]<-kle2cs_wt_wt_0mM - wt_wt_wt_0mM #kle2
  contrastsOI[["X.wt.wt.0mM_scc16cs_vs_wt"]]<-scc1cs_wt_wt_0mM - wt_wt_wt_0mM #scc1
  contrastsOI[["X.wt.wt.0mM_coh1cs_vs_wt"]]<-coh1cs_wt_wt_0mM - wt_wt_wt_0mM #coh1
  contrastsOI[["X.wt.wt.0mM_scc1coh1cs_vs_wt"]]<-scc1coh1cs_wt_wt_0mM - wt_wt_wt_0mM #scc1coh1
  contrastsOI[["X.wt.wt.0mM_scc1coh1cs-coh1"]]<-scc1coh1cs_wt_wt_0mM-coh1cs_wt_wt_0mM #scc1coh1cs-coh1
  contrastsOI[["X.wt.wt.0mM_scc1coh1cs-scc1"]]<-scc1coh1cs_wt_wt_0mM-scc1cs_wt_wt_0mM #scc1coh1cs-scc1

  contrastNames<-as.list(names(contrastsOI))
  names(contrastNames)<-c("aux","aux_sdc3BG","TIR1","TIR1aux","sdc3","dpy26","dpy26_sdc3BG","sdc3_dpy26BG","dpy26sdc3","aux_dpy26sdc3","kle2","scc1","coh1","scc1coh1","scc1coh1cs-coh1","scc1coh1cs-scc1")

  useContrasts<-c("aux_sdc3BG","sdc3","dpy26","dpy26sdc3","aux_dpy26sdc3","kle2","scc1","coh1","scc1coh1")
  DCsubset<-c("aux_sdc3BG","sdc3","dpy26","dpy26sdc3","aux_dpy26sdc3")
}


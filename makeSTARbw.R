library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library(GenomicFeatures)

source("functions.R")

###
### some variables
####
plotPDFs=F
fileNamePrefix="salmon_"
filterPrefix="noOsc_"
filterData=T
outPath="."
genomeVer="WS275"
genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)

genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6], IRanges(start=1, end=seqlengths(Celegans)[1:6]))

wbseqinfo<-seqinfo(Celegans)
seqnames(wbseqinfo)<-c(gsub("chr","",seqnames(Celegans)))
seqnames(wbseqinfo)<-c(gsub("^M$","MtDNA",seqnames(wbseqinfo)))
genome(wbseqinfo)<-genomeVer
ce11seqinfo<-seqinfo(Celegans)

makeDirs(outPath,dirNameList=c("tracks"))


fileList<-read.table(paste0(outPath,"/fastqList.txt"),stringsAsFactors=F,header=T)

# extract the strain variable
strain<-factor(as.character(unique(fileList$sampleName),levels=c("366","382","775","784")))
SMC<-strain
levels(SMC)<-c("wt","dpy26cs","kle2cs","scc1cs")

controlGrp<-levels(SMC)[1] # control group
groupsOI<-levels(SMC)[-1]


#######-
# Average STAR bigwig tracks ----------------------------------------------
######-

#biolSamples<-unique(sampleTable$SMC)
biolSamples<-c(grp,controlGrp)
avrFiles<-c()
for (biolSample in biolSamples) {
   idx<-which(sampleTable$SMC==biolSample)
   bwFiles<-paste0(outPath, "/tracks/",
                   sampleTable$sampleName[idx],
                   "_rpm.bw")
   wigFile<-paste0(outPath, "/tracks/", unique(sampleTable$SMC[idx]),
                   "_",unique(sampleTable$strain[idx]) ,"_rpm_Avr.wig")
   logwigFile<-gsub("_Avr","_logAvr", wigFile )

   finalOutputFile<-paste0(outPath,"/tracks/lfc_", grp, "_", controlGrp, "_ce11.bw")

   if(!file.exists(finalOutputFile)){
      system(paste0("wiggletools mean ", paste0(bwFiles,collapse=" "),
                    " > ", wigFile ))

      system(paste0("wiggletools offset 1 ", wigFile, " | wiggletools log 2 - >",
                    logwigFile))

      #wigToBigWig(x=wigFile, seqinfo=wbseqinfo,
      #             dest=gsub("\\.wig$","\\.bw",wigFile))
      #wigToBigWig(x=logwigFile, seqinfo=wbseqinfo,
      #            dest=gsub("\\.wig$","\\.bw",logwigFile))
      file.remove(wigFile)
   }
   avrFiles<-c(avrFiles,logwigFile)
}

# single log fold change track
if(!file.exists(finalOutputFile)){
   system(paste0("wiggletools diff ",paste0(avrFiles,collapse=" ")," > ",
                 paste0(outPath,"/tracks/lfc_",grp,"_",controlGrp,".wig")))
   wigToBigWig(x=paste0(outPath,"/tracks/lfc_",grp,"_",controlGrp,".wig"),
               seqinfo=wbseqinfo,
               dest=paste0(outPath,"/tracks/lfc_",grp,"_",controlGrp,".bw"))
   bw<-import(paste0(outPath,"/tracks/lfc_",grp,"_",controlGrp,".bw"))
   seqlevelsStyle(bw)<-"ucsc"
   idx<-match(seqlevels(ce11seqinfo),seqlevels(bw))
   seqinfo(bw,new2old=idx)<-ce11seqinfo
   export.bw(bw,finalOutputFile)
   file.remove(paste0(outPath,"/tracks/lfc_",grp,"_",controlGrp,".wig"))
   file.remove(paste0(outPath,"/tracks/lfc_",grp,"_",controlGrp,".bw"))
   file.remove(avrFiles)
}

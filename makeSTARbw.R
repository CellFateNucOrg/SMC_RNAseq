library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library(GenomicFeatures)

source("./functions.R")
source("./variableSettings.R")
###
### some variables
####

genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6], IRanges(start=1, end=seqlengths(Celegans)[1:6]))

wbseqinfo<-seqinfo(Celegans)
seqnames(wbseqinfo)<-c(gsub("chr","",seqnames(Celegans)))
seqnames(wbseqinfo)<-c(gsub("^M$","MtDNA",seqnames(wbseqinfo)))
genome(wbseqinfo)<-genomeVer
ce11seqinfo<-seqinfo(Celegans)

makeDirs(outPath,dirNameList=c("tracks"))




#######-
# Average STAR F&R bigwig tracks ----------------------------------------------
######-
if(remakeFiles){
   file.remove(list.files(path=paste0(outPath,"/tracks/"),pattern="STAR_lfc.*ce11\\.bw"))
}

finalFilesExist<-grep("STAR_lfc.*ce11\\.bw",list.files(paste0(outPath,"/tracks/")))

if(length(finalFilesExist)==0){
   for (i in 1:nrow(sampleTable)){
      # sum unique and multiunique forward
      wigFile<-paste0(outPath,"/tracks/STAR_",sampleTable$sampleName[i],"_F.wig")
      bwFiles<-paste0(outPath, "/bamSTAR/", sampleTable$sampleName[i],
                      c("_Unique_F.bw","_UniqueMultiple_F.bw"))
      system(paste0("wiggletools sum ", paste0(bwFiles,collapse=" "),
                    " > ", wigFile ))
      # sum unique and multiunique reverse
      wigFile<-paste0(outPath,"/tracks/STAR_",sampleTable$sampleName[i],"_R.wig")
      bwFiles<-paste0(outPath, "/bamSTAR/", sampleTable$sampleName[i],
                      c("_Unique_R.bw","_UniqueMultiple_R.bw"))
      system(paste0("wiggletools sum ", paste0(bwFiles,collapse=" "),
                    " > ", wigFile ))
      # average F & R
      wigFile<-paste0(outPath,"/tracks/STAR_",sampleTable$sampleName[i],".wig")
      wigFiles<-paste0(outPath,"/tracks/STAR_",sampleTable$sampleName[i],
                       c("_F.wig","_R.wig"))
      system(paste0("wiggletools mean ", paste0(wigFiles,collapse=" "),
                    " > ", wigFile ))
      file.remove(wigFiles)
   }

   #######-
   # Average STAR bigwig tracks by sample---------------------------------------
   ######-

   #biolSamples<-unique(sampleTable$SMC)
   for (grp in groupsOI) {
      biolSamples<-c(grp,controlGrp)
      avrFiles<-c()
      for (biolSample in biolSamples) {
         idx<-which(sampleTable$SMC==biolSample)
         bwFiles<-paste0(outPath, "/tracks/STAR_",
                         sampleTable$sampleName[idx],
                         ".wig")
         wigFile<-paste0(outPath, "/tracks/STAR_", unique(sampleTable$SMC[idx]),
                         "_",unique(sampleTable$strain[idx]) ,"_Avr.wig")
         logwigFile<-gsub("_Avr","_logAvr", wigFile )

         finalOutputFile<-paste0(outPath,"/tracks/STAR_lfc_", grp, "_", controlGrp, "_ce11.bw")

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
            file.remove(bwFiles)
         }
         avrFiles<-c(avrFiles,logwigFile)
      }

      # single log fold change track
      if(!file.exists(finalOutputFile)){
         system(paste0("wiggletools diff ",paste0(avrFiles,collapse=" ")," > ",
                       paste0(outPath,"/tracks/lfc_",grp,"_",controlGrp,".wig")))
         rtracklayer::wigToBigWig(x=paste0(outPath,"/tracks/lfc_",grp,"_",controlGrp,".wig"),
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
   }
}

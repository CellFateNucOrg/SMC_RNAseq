library("rjson")
library(dplyr)
## get number of reads at each stage of process
source("./variableSettings.R")

if(!file.exists(paste0(outPath,"/qc/rawData/readCount.txt"))){
  system(paste0("./collectFastQCreadCounts.sh ",outPath,"/qc/rawData"))
}

if(!file.exists(paste0(outPath,"/qc/cutadapt/readCount.txt"))){
  system(paste0("./collectFastQCreadCounts.sh ",outPath,"/qc/cutadapt"))
}

rawDataCount<-read.delim(paste0(outPath,"/qc/rawData/readCount.txt"),header=F)
colnames(rawDataCount)<-c("step","library","fragCount")
cutadaptCount<-read.delim(paste0(outPath,"/qc/cutadapt/readCount.txt"),header=F)
colnames(cutadaptCount)<-c("step","library","fragCount")
rawDataCount$fragCount==cutadaptCount$fragCount

countTable<-data.frame(library=gsub("_fastqc$","",cutadaptCount$library), rawData=rawDataCount$fragCount,
                       cutadapt=cutadaptCount$fragCount,
                       salmonNumMap=NA, salmonPercentMap=NA,
                       starUniqMap=NA, starPercentUniqMap=NA,
                       MreadsPerLib=NA)


salmonDirs<-list.files(path=paste0(outPath,"/salmon/mRNA"))
d=salmonDirs[1]
for (d in salmonDirs){
  jsonData<-fromJSON(file=paste0(outPath,"/salmon/mRNA/",d,"/aux_info/meta_info.json"))
  i<-which(countTable$library==d)
  countTable$salmonNumMap[i]<-jsonData$num_mapped
  countTable$salmonPercentMap[i]<-round(jsonData$percent_mapped,2)
}


starFiles<-list.files(path=paste0(outPath,"/bamSTAR"),pattern="_Log\\.final\\.out$")

for (f in starFiles) {
  df<-read.delim(paste0(outPath,"/bamSTAR/",f),header=F,stringsAsFactors=F)
  i<-grep("Uniquely mapped reads number",df$V1)
  j<-grep("Uniquely mapped reads %",df$V1)
  lib<-gsub("_Log\\.final\\.out","",f)
  countTable[countTable$library==lib,"starUniqMap"]<-df$V2[i]
  countTable[countTable$library==lib,"starPercentUniqMap"]<-df$V2[j]
}


#countTable$percentRawMapped<-round(100*countTable$salmonNumMap/countTable$rawData,2)

#seqlib<-gsub("_L[1|2]$","",countTable$library)


write.table(countTable,file=paste0(outPath,"/qc/readCountsByStage.txt"), col.names=T, row.names=F, quote=F)

countTable$seqlib<-as.factor(gsub("_L[1|2]$","",countTable$library))
by_lib<-countTable %>% group_by(seqlib) %>% summarise(total=round(sum(rawData)/1e6,1))
min(by_lib$total)
max(by_lib$total)

library("rjson")

## get number of reads at each stage of process
outPath="."

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
			salmonNumMap_ncRNA=NA, salmonPercentMap_ncRNA=NA,
			salmonNumMap_tnRNA=NA, salmonPercentMap_tnRNA=NA,
			salmonNumMap_pseudoRNA=NA, salmonPercentMap_pseudoRNA=NA,
			salmonNumMap_rptRNA=NA, salmonPercentMap_rptRNA=NA,
                       starUniqMap=NA, starPercentUniqMap=NA,
			starUniqMap_rpt=NA, starPercentUniqMap_rpt=NA,
			bwaUniqMap_random=NA,	bwaPercentUniqMap_random=NA,
			bwaUniqMap_none=NA,   bwaPercentUniqMap_none=NA)


salmonDirs<-list.files(path=paste0(outPath,"/salmon/mRNA"))
d=salmonDirs[1]
for (d in salmonDirs){
  jsonData<-fromJSON(file=paste0(outPath,"/salmon/mRNA/",d,"/aux_info/meta_info.json"))
  i<-which(countTable$library==d)
  countTable$salmonNumMap[i]<-jsonData$num_mapped
  countTable$salmonPercentMap[i]<-round(jsonData$percent_mapped,2)
}

salmonDirs<-list.files(path=paste0(outPath,"/salmon/ncRNA"))
d=salmonDirs[1]
for (d in salmonDirs){
  jsonData<-fromJSON(file=paste0(outPath,"/salmon/ncRNA/",d,"/aux_info/meta_info.json"))
  i<-which(countTable$library==d)
  countTable$salmonNumMap_ncRNA[i]<-jsonData$num_mapped
  countTable$salmonPercentMap_ncRNA[i]<-round(jsonData$percent_mapped,2)
}

salmonDirs<-list.files(path=paste0(outPath,"/salmon/tnRNA"))
d=salmonDirs[1]
for (d in salmonDirs){
  jsonData<-fromJSON(file=paste0(outPath,"/salmon/tnRNA/",d,"/aux_info/meta_info.json"))
  i<-which(countTable$library==d)
  countTable$salmonNumMap_tnRNA[i]<-jsonData$num_mapped
  countTable$salmonPercentMap_tnRNA[i]<-round(jsonData$percent_mapped,2)
}

salmonDirs<-list.files(path=paste0(outPath,"/salmon/pseudoRNA"))
d=salmonDirs[1]
for (d in salmonDirs){
  jsonData<-fromJSON(file=paste0(outPath,"/salmon/pseudoRNA/",d,"/aux_info/meta_info.json"))
  i<-which(countTable$library==d)
  countTable$salmonNumMap_pseudoRNA[i]<-jsonData$num_mapped
  countTable$salmonPercentMap_pseudoRNA[i]<-round(jsonData$percent_mapped,2)
}

salmonDirs<-list.files(path=paste0(outPath,"/salmon/rptRNA_15"))
d=salmonDirs[1]
for (d in salmonDirs){
  jsonData<-fromJSON(file=paste0(outPath,"/salmon/rptRNA_15/",d,"/aux_info/meta_info.json"))
  i<-which(countTable$library==d)
  countTable$salmonNumMap_rptRNA[i]<-jsonData$num_mapped
  countTable$salmonPercentMap_rptRNA[i]<-round(jsonData$percent_mapped,2)
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


starFiles<-list.files(path=paste0(outPath,"/bamSTARrpt"),pattern="_Log\\.final\\.out$")

for (f in starFiles) {
  df<-read.delim(paste0(outPath,"/bamSTARrpt/",f),header=F,stringsAsFactors=F)
  i<-grep("Uniquely mapped reads number",df$V1)
  j<-grep("Uniquely mapped reads %",df$V1)
  lib<-gsub("_Log\\.final\\.out","",f)
  countTable[countTable$library==lib,"starUniqMap_rpt"]<-df$V2[i]
  countTable[countTable$library==lib,"starPercentUniqMap_rpt"]<-df$V2[j]
}



bwaFiles<-list.files(path=paste0(outPath,"/bamBWA"),pattern="_union_random\\.txt$")
for (f in bamFiles) {
  df<-read.delim(paste0(outPath,"/bamBWA/",f),header=F,stringsAsFactors=F)
  lib<-gsub("_union_random\\.txt","",f)
  countTable[countTable$library==lib,"bwaUniqMap_random"]<-sum(df$V2)
}
countTable$bwaPercentUniqMap_random=100*countTable$bwaUniqMap_random/countTable$cutadapt


bwaFiles<-list.files(path=paste0(outPath,"/bamBWA"),pattern="_union_none\\.txt$")
for (f in bamFiles) {
  df<-read.delim(paste0(outPath,"/bamBWA/",f),header=F,stringsAsFactors=F)
  lib<-gsub("_union_none\\.txt","",f)
  countTable[countTable$library==lib,"bwaUniqMap_none"]<-sum(df$V2)
}

countTable$bwaPercentUniqMap_none=100*countTable$bwaUniqMap_none/countTable$cutadapt


#countTable$percentRawMapped<-round(100*countTable$salmonNumMap/countTable$rawData,2)

write.table(countTable,file=paste0(outPath,"/qc/collectedReadCounts.txt"), col.names=T, row.names=F, quote=F)

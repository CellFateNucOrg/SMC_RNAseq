library(readxl)
library(ggVennDiagram)
library(ggplot2)
library(EnhancedVolcano)
library(wormcat)
library(xlsx)

source("functions.R")
source("./variableSettings.R")

fileList<-read.table(paste0(outPath,"/fastqList.txt"),stringsAsFactors=F,header=T)

# extract the strain variable
strain<-factor(as.character(unique(fileList$sampleName),levels=c("366","382","775","784")))
SMC<-strain
levels(SMC)<-c("wt","dpy26cs","kle2cs","scc1cs")

controlGrp<-levels(SMC)[1] # control group
groupsOI<-levels(SMC)[-1]


############################
## Gene name lists for wormcat
############################
## need to process by dropping excel file into http://www.wormcat.com
## since r function does not seem to work.

if(!dir.exists(paste0(outPath,"/wormcat/p",padjVal,"_lfc",lfcVal,"/"))) {
  dir.create(paste0(outPath,"/wormcat/p",padjVal,"_lfc",lfcVal,"/"))
}
# ## all genes
# sigTables<-list()
# for (grp in groupsOI){
#   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))
#
#   sigTables[[paste0(grp,"_all")]]<-as.data.frame(
#     getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
#                           namePadjCol="padj",
#                           nameLfcCol="log2FoldChange",
#                           direction="both",
#                           chr="all", nameChrCol="chr"))
# }
# sigGenes<-lapply(sigTables, "[", ,"wormbaseID")




### upregulated genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

  sigTables[[paste0(grp,"_up")]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="gt",
                        chr="all", nameChrCol="chr"))
}
sigGenesUp<-lapply(sigTables, "[", ,"wormbaseID")



### down regulated genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

  sigTables[[paste0(grp,"_down")]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=-lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="lt",
                        chr="all", nameChrCol="chr"))

}

sigGenesDown<-lapply(sigTables, "[", ,"wormbaseID")


wormcatIn<-c(sigGenesUp,sigGenesDown)

wormcatIn<-lapply(wormcatIn,function (x) c("Wormbase ID",x))

openxlsx::write.xlsx(wormcatIn,
                     file=paste0(outPath,"/wormcat/",fileNamePrefix,"wormcat.xlsx"))


#worm_cat_fun( file_to_process=paste0(outPath,"/wormcat/wormcat.xlsx"), output_dir=paste0(outPath,"/wormcat/"), input_type="Wormbase.ID")

#if(!dir.exists(paste0(outPath,"/tissue")) { dir.create(paste0(outPath,"/tissue")) }

#####################
## worm tissue
#####################
#http://worm-tissue.princeton.edu/search
# need lists of sequence IDs.

# read in tissue prediction data set to screen out entrez ids not included
if(!file.exists("all_tissue_prediction_scores.txt")){
  system("wget https://worm.princeton.edu/media/download/all_tissue_prediction_scores.txt")
}
tissueScores<-read.delim("./all_tissue_prediction_scores.txt")
write.table(tissueScores$entrez,"tissueScoresEntrezIDs.txt",quote=F,
            row.names=F,col.names=F)

problemIDs<-c(177376, 190779, 188651, 13217889, 3565958, 13183957, 189501, 185871, 6418805, 178591, 183939, 172176, 185960, 186971, 177183, 185400, 187634, 181346, 246001, 179099, 178276, 189582, 4926949, 13192684, 173066, 176017, 13218347, 189137, 13181785, 266818, 184170, 191361, 3896793, 13186920, 183376, 13180879, 186190, 181938, 191092, 6418577, 190412, 176677, 188272, 190032, 3565863, 174389, 190519, 183030, 260149, 3565471, 171717, 3896778, 13213331, 13218463, 7040166, 13191473, 13200771, 190944, 6418869, 13216772, 3896811, 185914, 191226, 3565324, 13208938, 13187086, 353398, 183874, 190152, 188264, 13183020, 172113, 172069, 184874, 3565662, 353484, 3565976, 190548, 190881, 191193, 13203945, 13198804, 183620, 182872, 189678, 190593, 188822, 181049, 185401, 3565900, 183099, 184649, 13192683, 190322, 4363122, 173636, 184894, 183302, 175224, 3896751, 179697, 13182022, 182928, 13181792, 184592, 177322, 178497, 13184661, 176248, 177062, 178705, 183727, 191080, 186316, 172271, 182522, 190652, 3564771, 3896762, 190591, 185363, 13214333, 184415, 353379, 184544, 187752, 3565988, 3565306, 260236, 191531, 3565755, 3565094, 13186136, 182968, 13211145, 3565657, 13189669, 353432, 178539, 3565309, 188411, 190535, 266943, 3564963, 190265, 183920, 13198671, 185398, 13181788, 188721, 3565350, 13218303, 188670, 178610, 190517, 173570, 186124, 190599, 3896890, 24104801, 184390, 181419, 187319, 3565173, 177127, 185795, 182828, 13198805, 188220, 13182910, 182034, 179317, 191530, 190042, 191289, 191512, 188727, 3896846, 189550, 181709, 175746, 189344, 190785, 188640, 189775, 171600, 191082, 3565741, 179799, 3565119, 3565523, 173318, 190191, 182583, 259553, 13183522, 176535, 3565948, 3565913, 180903, 13182666, 183459, 183872, 3565023, 190924, 266910, 190318, 13185431, 185710, 13189099, 24104176, 185420, 188516, 3565870, 189105, 3565963, 3565412, 3564959, 13218911, 189136, 189549, 183852, 189493, 186733, 3565844, 3565671, 172043, 189133, 176529, 4926974, 180471, 183050, 183058, 188636, 190628, 183518, 182282, 190654, 13179148, 177239, 13180978, 184598, 188263, 189774, 3896770, 179982, 6418806, 13217482, 190034, 177025, 188517, 353435, 178608, 187701, 187896, 189351, 177156, 185943, 3896853, 3565319, 13224779, 181988, 184593, 353485, 182510, 13183021, 13179162, 175127, 172272)



if(!dir.exists(paste0(outPath,"/tissue/wormtissue/p",padjVal,"_lfc",lfcVal,"/"))) {
  dir.create(paste0(outPath,"/tissue/wormtissue/p",padjVal,"_lfc",lfcVal,"/"))
}

# ## significantly changed genes
# sigTables<-list()
# for (grp in groupsOI){
#   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))
#
#   sigTables[[paste0(grp)]]<-as.data.frame(
#     getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
#                         namePadjCol="padj",
#                         nameLfcCol="log2FoldChange",
#                         direction="both",
#                         chr="all", nameChrCol="chr"))
# }
# sigGenes<-lapply(sigTables, "[", ,"entrezID")
# lapply(sigGenes,length)
# sigGenes<-lapply(sigGenes,na.omit)
#
# for (grp in groupsOI){
#   subset<-sigGenes[[grp]][sigGenes[[grp]] %in% tissueScores$entrez & !( sigGenes[[grp]] %in% problemIDs) ]
#   print(paste(grp,length(subset),"genes"))
#   write.table(subset, file=paste0(outPath,"/tissue/wormtissue/",grp,"_allGenes_ENTREZ.txt"), quote=F, row.names=F,col.names=F)
# }


## upregulated genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

  sigTables[[paste0(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="gt",
                        chr="all", nameChrCol="chr"))
}
sigGenes<-lapply(sigTables, "[", ,"entrezID")
lapply(sigGenes,length)
sigGenes<-lapply(sigGenes,na.omit)

for (grp in groupsOI){
  subset<-sigGenes[[grp]][sigGenes[[grp]] %in% tissueScores$entrez & !( sigGenes[[grp]] %in% problemIDs) ]
  print(paste(grp,length(subset),"genes"))
  write.table(subset, file=paste0(outPath,"/tissue/wormtissue/",fileNamePrefix,grp,"_upGenes_ENTREZ.txt"), quote=F, row.names=F,col.names=F)
}


## downregulated genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

  sigTables[[paste0(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=-lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="lt",
                        chr="all", nameChrCol="chr"))
}
sigGenes<-lapply(sigTables, "[", ,"entrezID")
lapply(sigGenes,length)
sigGenes<-lapply(sigGenes,na.omit)

for (grp in groupsOI){
  subset<-sigGenes[[grp]][sigGenes[[grp]] %in% tissueScores$entrez & !( sigGenes[[grp]] %in% problemIDs) ]
  print(paste(grp,length(subset),"genes"))
  write.table(subset, file=paste0(outPath,"/tissue/wormtissue/",fileNamePrefix,grp,"_downGenes_ENTREZ.txt"), quote=F, row.names=F,col.names=F)
}



####### sequenceID

## upregulated genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

  sigTables[[paste0(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="gt",
                        chr="all", nameChrCol="chr"))
}
sigGenes<-lapply(sigTables, "[", ,"sequenceID")
lapply(sigGenes,length)
sigGenes<-lapply(sigGenes,na.omit)

for (grp in groupsOI){
  print(paste(grp,length(sigGenes[[grp]]),"genes"))
  write.table(sigGenes[[grp]], file=paste0(outPath,"/tissue/wormtissue/",fileNamePrefix,grp,"_upGenes_sequenceID.txt"), quote=F, row.names=F,col.names=F)
}


## downregulated genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

  sigTables[[paste0(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=-lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="lt",
                        chr="all", nameChrCol="chr"))
}
sigGenes<-lapply(sigTables, "[", ,"sequenceID")
lapply(sigGenes,length)
sigGenes<-lapply(sigGenes,na.omit)

for (grp in groupsOI){
  print(paste(grp,length(sigGenes[[grp]]),"genes"))
  write.table(sigGenes[[grp]], file=paste0(outPath,"/tissue/wormtissue/",fileNamePrefix,grp,"_downGenes_sequenceID.txt"), quote=F, row.names=F,col.names=F)
}



########
## TEA - tissue enrichment analysis
########
# https://www.micropublication.org/media/2018/03/microPublication.biology-10.17912-W25Q2N.pdf

if(!dir.exists(paste0(outPath,"/tissue/tea/p",padjVal,"_lfc",lfcVal,"/"))) {
  dir.create(paste0(outPath,"/tissue/tea/p",padjVal,"_lfc",lfcVal,"/"))
}


# ## significantly changed genes
# sigTables<-list()
# for (grp in groupsOI){
#   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))
#
#   sigTables[[paste0(grp)]]<-as.data.frame(
#     getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
#                         namePadjCol="padj",
#                         nameLfcCol="log2FoldChange",
#                         direction="both",
#                         chr="all", nameChrCol="chr"))
# }
# sigGenes<-lapply(sigTables, "[", ,"wormbaseID")
# lapply(sigGenes,length)
# sigGenes<-lapply(sigGenes,na.omit)


sink(file=paste0(outPath,"/runTea.sh"),append=FALSE, type="output")
cat("#! /bin/bash\n")
cat(paste0("cd ./tissue/tea/p",padjVal,"_lfc",lfcVal,"\n"))
#cat("cd ./tissue/tea\n")
sink()
# #file.create(paste0(outPath,"/runTea.sh"),overwrite=T)
# for (grp in groupsOI){
#   write.table(sigGenes[[grp]], file=paste0(outPath,"/tissue/tea/",grp,"_allGenes_WBID.txt"), quote=F, row.names=F,col.names=F)
#   sink(file=paste0(outPath,"/runTea.sh"),append=TRUE, type="output")
#   cat(paste0("tea -q 0.05 -s ",grp,"_allGenes_WBID.txt ", grp,"_all_tissue tissue\n"))
#   cat(paste0("tea -q 0.05 -s ",grp,"_allGenes_WBID.txt ", grp,"_all_phe phenotype\n"))
#   cat(paste0("tea -q 0.05 -s ",grp,"_allGenes_WBID.txt ", grp,"_all_go go\n"))
#   sink()
# }

## upregulated genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

  sigTables[[paste0(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="gt",
                        chr="all", nameChrCol="chr"))
}
sigGenes<-lapply(sigTables, "[", ,"wormbaseID")
lapply(sigGenes,length)
sigGenes<-lapply(sigGenes,na.omit)

partialPrefix=gsub(paste0("p",padjVal,"_lfc",lfcVal,"/"),"",fileNamePrefix)
for (grp in groupsOI){
  write.table(sigGenes[[grp]], file=paste0(outPath,"/tissue/tea/",fileNamePrefix,grp,"_upGenes_WBID.txt"), quote=F, row.names=F,col.names=F)
  sink(file=paste0(outPath,"/runTea.sh"),append=TRUE, type="output")
  cat(paste0("tea -q 0.05 -s ",partialPrefix,grp,"_upGenes_WBID.txt ", partialPrefix,grp,"_up_tissue tissue\n"))
  cat(paste0("tea -q 0.05 -s ",partialPrefix,grp,"_upGenes_WBID.txt ", partialPrefix,grp,"_up_phe phenotype\n"))
  cat(paste0("tea -q 0.05 -s ",partialPrefix,grp,"_upGenes_WBID.txt ", partialPrefix,grp,"_up_go go\n"))
  sink()
}


## downregulated genes
sigTables<-list()
for (grp in groupsOI){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,grp,"_DESeq2_fullResults.rds"))

  sigTables[[paste0(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=-lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="lt",
                        chr="all", nameChrCol="chr"))
}
sigGenes<-lapply(sigTables, "[", ,"wormbaseID")
lapply(sigGenes,length)
sigGenes<-lapply(sigGenes,na.omit)

for (grp in groupsOI){
  write.table(sigGenes[[grp]], file=paste0(outPath,"/tissue/tea/",fileNamePrefix,grp,"_downGenes_WBID.txt"), quote=F, row.names=F,col.names=F)
  sink(file=paste0(outPath,"/runTea.sh"),append=TRUE, type="output")
  cat(paste0("tea -q 0.05 -s ",partialPrefix,grp,"_downGenes_WBID.txt ", partialPrefix, grp,"_down_tissue tissue\n"))
  cat(paste0("tea -q 0.05 -s ",partialPrefix,grp,"_downGenes_WBID.txt ", partialPrefix ,grp,"_down_phe phenotype\n"))
  cat(paste0("tea -q 0.05 -s ",partialPrefix,grp,"_downGenes_WBID.txt ", partialPrefix, grp, "_down_go go\n"))
  sink()
}

sink(file=paste0(outPath,"/runTea.sh"),append=TRUE, type="output")
cat("cd ../../../\n")
sink()

system(paste0("chmod +x ",outPath,"/runTea.sh"))
system(paste0(outPath,"/runTea.sh")) # doesn't work? but can be run from command line

library(readxl)
library(ggplot2)
library(wormcat)
library(xlsx)
library(reticulate)
library(dplyr)
#conda_install(envname="tea",packages="tissue_enrichment_analysis",pip=T, pip_options="git+https://github.com/dangeles/TissueEnrichmentAnalysis.git")
## note: there is a bug in this version of tea which need to be corrected in the
## code. Go to this file: vi ~/miniconda3/envs/tea/bin/tea and on line 66 change args.tisse_dictionary to args.dictionary
# install wormcat from my fork on github devtools::install_github("jsemple19/Wormcat")

source("functions.R")
source("./variableSettings.R")
if(filterData){
  fileNamePrefix<-filterPrefix
}
metadata<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))


############################-
## Wormcat-----
############################-
## need to process by dropping excel file into http://www.wormcat.com
## since r function does not seem to work.

if(!dir.exists(paste0(outPath,"/wormcat/p",padjVal,"_lfc",lfcVal,"/"))) {
  dir.create(paste0(outPath,"/wormcat/p",padjVal,"_lfc",lfcVal,"/"),
             recursive=T)
}

wormcatDataURL="http://www.wormcat.com/static/download/whole_genome_nov-16-2019.csv"
wormcatData=paste0(getwd(),"/publicData/",basename(wormcatDataURL))
if(!file.exists(wormcatData)){
  download.file(url=wormcatDataURL,destfile=wormcatData)
}
#read.csv(wormcatData)
# ## all genes
# sigTables<-list()
# for (grp in useContrasts){
#   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
#
#   sigTables[[paste0(grp,"_all")]]<-as.data.frame(
#     getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
#                           namePadjCol="padj",
#                           nameLfcCol="log2FoldChange",
#                           direction="both",
#                           chr="all", nameChrCol="chr"))
# }
# sigGenes<-lapply(sigTables, "[", ,"wormbaseID")



for(grp in useContrasts){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

  ### upregulated genes
  sigTable<-data.frame(Wormbase.ID=getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="gt",
                        chr="all", nameChrCol="chr")$wormbaseID)
  write.table(sigTable,paste0(outPath,"/wormcat/",fileNamePrefix,contrastNames[[grp]],
                              "_up_wormcat.csv"),
            row.names=F,quote=F,col.names=T)
  worm_cat_fun(file_to_process=paste0(outPath,"/wormcat/",fileNamePrefix,contrastNames[[grp]],
                                       "_up_wormcat.csv"),
                title=paste(grp,"up"),
                output_dir=paste0(outPath,"/wormcat/",fileNamePrefix,contrastNames[[grp]],"_up"),
                annotation_file=wormcatData,
                input_type="Wormbase.ID", rm_dir=FALSE,
                zip_files=FALSE)

  ### down regulated genes
  sigTable<-data.frame(Wormbase.ID=getSignificantGenes(salmon, padj=padjVal, lfc= -lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="lt",
                        chr="all", nameChrCol="chr")$wormbaseID)
  write.table(sigTable,paste0(outPath,"/wormcat/",fileNamePrefix,contrastNames[[grp]],"_down_wormcat.csv"),
              row.names=F,quote=F,col.names=T)
  worm_cat_fun( file_to_process=paste0(outPath,"/wormcat/", fileNamePrefix, contrastNames[[grp]],
                                       "_down_wormcat.csv"),
                title=paste(grp,"down"),
                output_dir=paste0(outPath,"/wormcat/",fileNamePrefix,contrastNames[[grp]],"_down"),
                annotation_file=wormcatData,
                input_type="Wormbase.ID", rm_dir=FALSE,
                zip_files=FALSE)
}






#if(!dir.exists(paste0(outPath,"/tissue")) { dir.create(paste0(outPath,"/tissue")) }

#####################-
## worm tissue-----
#####################-
#http://worm-tissue.princeton.edu/search
#http://worm.princeton.edu/
# need lists of sequence IDs.

# read in tissue prediction data set to screen out entrez ids not included
tissuePredScores<-paste0(outPath,"/publicData/all_tissue_prediction_scores.txt")
if(remakeFiles & file.exists(tissuePredScores)){
  file.remove(tissuePredScores)
}
if(!file.exists(tissuePredScores)){
  download.file(url="https://worm.princeton.edu/media/download/all_tissue_prediction_scores.txt", destfile=tissuePredScores)
}
tissueScores<-read.delim(tissuePredScores)
write.table(tissueScores$entrez,
            paste0(outPath,"/publicData/tissueScoresEntrezIDs.txt"),
            quote=F,row.names=F,col.names=F)

problemIDs<-c(177376, 190779, 188651, 13217889, 3565958, 13183957, 189501, 185871, 6418805, 178591, 183939, 172176, 185960, 186971, 177183, 185400, 187634, 181346, 246001, 179099, 178276, 189582, 4926949, 13192684, 173066, 176017, 13218347, 189137, 13181785, 266818, 184170, 191361, 3896793, 13186920, 183376, 13180879, 186190, 181938, 191092, 6418577, 190412, 176677, 188272, 190032, 3565863, 174389, 190519, 183030, 260149, 3565471, 171717, 3896778, 13213331, 13218463, 7040166, 13191473, 13200771, 190944, 6418869, 13216772, 3896811, 185914, 191226, 3565324, 13208938, 13187086, 353398, 183874, 190152, 188264, 13183020, 172113, 172069, 184874, 3565662, 353484, 3565976, 190548, 190881, 191193, 13203945, 13198804, 183620, 182872, 189678, 190593, 188822, 181049, 185401, 3565900, 183099, 184649, 13192683, 190322, 4363122, 173636, 184894, 183302, 175224, 3896751, 179697, 13182022, 182928, 13181792, 184592, 177322, 178497, 13184661, 176248, 177062, 178705, 183727, 191080, 186316, 172271, 182522, 190652, 3564771, 3896762, 190591, 185363, 13214333, 184415, 353379, 184544, 187752, 3565988, 3565306, 260236, 191531, 3565755, 3565094, 13186136, 182968, 13211145, 3565657, 13189669, 353432, 178539, 3565309, 188411, 190535, 266943, 3564963, 190265, 183920, 13198671, 185398, 13181788, 188721, 3565350, 13218303, 188670, 178610, 190517, 173570, 186124, 190599, 3896890, 24104801, 184390, 181419, 187319, 3565173, 177127, 185795, 182828, 13198805, 188220, 13182910, 182034, 179317, 191530, 190042, 191289, 191512, 188727, 3896846, 189550, 181709, 175746, 189344, 190785, 188640, 189775, 171600, 191082, 3565741, 179799, 3565119, 3565523, 173318, 190191, 182583, 259553, 13183522, 176535, 3565948, 3565913, 180903, 13182666, 183459, 183872, 3565023, 190924, 266910, 190318, 13185431, 185710, 13189099, 24104176, 185420, 188516, 3565870, 189105, 3565963, 3565412, 3564959, 13218911, 189136, 189549, 183852, 189493, 186733, 3565844, 3565671, 172043, 189133, 176529, 4926974, 180471, 183050, 183058, 188636, 190628, 183518, 182282, 190654, 13179148, 177239, 13180978, 184598, 188263, 189774, 3896770, 179982, 6418806, 13217482, 190034, 177025, 188517, 353435, 178608, 187701, 187896, 189351, 177156, 185943, 3896853, 3565319, 13224779, 181988, 184593, 353485, 182510, 13183021, 13179162, 175127, 172272)



if(!dir.exists(paste0(outPath,"/tissue/wormtissue/p",padjVal,"_lfc",lfcVal,"/"))) {
  dir.create(paste0(outPath,"/tissue/wormtissue/p",padjVal,"_lfc",lfcVal,"/"),
             recursive=T)
}


## upregulated genes
sigTables<-list()
for (grp in useContrasts){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  salmon<-salmon[salmon$entrezID %in% tissueScores$entrez,]
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

for(grp in useContrasts){
  subset<-sigGenes[[grp]][sigGenes[[grp]] %in% tissueScores$entrez & !( sigGenes[[grp]] %in% problemIDs) ]
  print(paste(grp,length(subset),"genes"))
  write.table(subset, file=paste0(outPath,"/tissue/wormtissue/",fileNamePrefix,contrastNames[[grp]],"_upGenes_ENTREZ.txt"), quote=F, row.names=F,col.names=F)
}


## downregulated genes
sigTables<-list()
for(grp in useContrasts){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  salmon<-salmon[salmon$entrezID %in% tissueScores$entrez,]
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

for(grp in useContrasts){
  subset<-sigGenes[[grp]][sigGenes[[grp]] %in% tissueScores$entrez & !( sigGenes[[grp]] %in% problemIDs) ]
  print(paste(grp,length(subset),"genes"))
  write.table(subset, file=paste0(outPath,"/tissue/wormtissue/",fileNamePrefix,contrastNames[[grp]],"_downGenes_ENTREZ.txt"), quote=F, row.names=F,col.names=F)
}



####### sequenceID-----

## upregulated genes
sigTables<-list()
for(grp in useContrasts){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  salmon<-salmon[salmon$entrezID %in% tissueScores$entrez,]
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

for(grp in useContrasts){
  print(paste(grp,length(sigGenes[[grp]]),"genes"))
  write.table(sigGenes[[grp]], file=paste0(outPath,"/tissue/wormtissue/",
                        fileNamePrefix,contrastNames[[grp]],
                        "_upGenes_sequenceID.txt"),
              quote=F, row.names=F,col.names=F)
}


## downregulated genes
sigTables<-list()
for(grp in useContrasts){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
  salmon<-salmon[salmon$entrezID %in% tissueScores$entrez,]
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

for(grp in useContrasts){
  print(paste(grp,length(sigGenes[[grp]]),"genes"))
  write.table(sigGenes[[grp]], file=paste0(outPath,"/tissue/wormtissue/",fileNamePrefix,contrastNames[[grp]],"_downGenes_sequenceID.txt"), quote=F, row.names=F,col.names=F)
}



########-
## TEA - tissue enrichment analysis-----
########-
# https://www.micropublication.org/media/2018/03/microPublication.biology-10.17912-W25Q2N.pdf

if(!dir.exists(paste0(outPath,"/tissue/tea/p",padjVal,"_lfc",lfcVal,"/"))) {
  dir.create(paste0(outPath,"/tissue/tea/p",padjVal,"_lfc",lfcVal,"/"),
             recursive=T)
}


# fetch dictionaries
# anatomy
anaURL=paste0("http://caltech.wormbase.org/TissueEnrichmentAnalysis/DICTs/anatomy_dict_95_33_", genomeVer, ".csv")
anaDict=paste0(outPath, "/publicData/",basename(anaURL))
if(!file.exists(anaDict)){
  download.file(url=anaURL, destfile=anaDict)
}

# go
goURL=paste0("http://caltech.wormbase.org/TissueEnrichmentAnalysis/DICTs/go_dict_95_100_", genomeVer, ".csv")
goDict=paste0(outPath, "/publicData/",basename(goURL))
if(!file.exists(goDict)){
  download.file(url=goURL, destfile=goDict)
}

# phenotype
pheURL=paste0("http://caltech.wormbase.org/TissueEnrichmentAnalysis/DICTs/phenotype_dict_95_50_", genomeVer, ".csv")
pheDict=paste0(outPath,"/publicData/",basename(pheURL))
if(!file.exists(pheDict)){
  download.file(url=pheURL, destfile=pheDict)
}
# ## significantly changed genes
# sigTables<-list()
# for (grp in useContrasts){
#   salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))
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


condaActivate<-gsub("conda$","activate",conda_binary(conda = "auto"))
sink(file=paste0(outPath,"/runTea.sh"),append=FALSE, type="output")
cat("#! /bin/bash\n")
cat(paste0("source ",condaActivate, " tea\n"))
cat(paste0("cd ./tissue/tea/p",padjVal,"_lfc",lfcVal,"\n"))
#cat("cd ./tissue/tea\n")
sink()
# #file.create(paste0(outPath,"/runTea.sh"),overwrite=T)
# for (grp in useContrasts){
#   write.table(sigGenes[[grp]], file=paste0(outPath,"/tissue/tea/",contrastNames[[grp]],"_allGenes_WBID.txt"), quote=F, row.names=F,col.names=F)
#   sink(file=paste0(outPath,"/runTea.sh"),append=TRUE, type="output")
#   cat(paste0("tea -q 0.05 -s ",contrastNames[[grp]],"_allGenes_WBID.txt ", contrastNames[[grp]],"_all_tissue tissue\n"))
#   cat(paste0("tea -q 0.05 -s ",contrastNames[[grp]],"_allGenes_WBID.txt ", contrastNames[[grp]],"_all_phe phenotype\n"))
#   cat(paste0("tea -q 0.05 -s ",contrastNames[[grp]],"_allGenes_WBID.txt ", contrastNames[[grp]],"_all_go go\n"))
#   sink()
# }

## upregulated genes
sigTables<-list()
for(grp in useContrasts){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

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
for(grp in useContrasts){
  write.table(sigGenes[[grp]], file=paste0(outPath,"/tissue/tea/",fileNamePrefix,contrastNames[[grp]],"_upGenes_WBID.txt"), quote=F, row.names=F,col.names=F)
  sink(file=paste0(outPath,"/runTea.sh"),append=TRUE, type="output")
  cat(paste0("tea -d ../../../",anaDict," -q 0.05 -s ",partialPrefix,contrastNames[[grp]],"_upGenes_WBID.txt ", partialPrefix,contrastNames[[grp]],"_up_tissue tissue\n"))
  cat(paste0("tea -d ../../../",pheDict," -q 0.05 -s ",partialPrefix,contrastNames[[grp]],"_upGenes_WBID.txt ", partialPrefix,contrastNames[[grp]],"_up_phe phenotype\n"))
  cat(paste0("tea -d ../../../",goDict," -q 0.05 -s ",partialPrefix,contrastNames[[grp]],"_upGenes_WBID.txt ", partialPrefix,contrastNames[[grp]],"_up_go go\n"))
  sink()
}


## downregulated genes
sigTables<-list()
for(grp in useContrasts){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

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

for(grp in useContrasts){
  write.table(sigGenes[[grp]], file=paste0(outPath,"/tissue/tea/",fileNamePrefix,contrastNames[[grp]],"_downGenes_WBID.txt"), quote=F, row.names=F,col.names=F)
  sink(file=paste0(outPath,"/runTea.sh"),append=TRUE, type="output")
  cat(paste0("tea -d ../../../",anaDict," -q 0.05 -s ",partialPrefix,contrastNames[[grp]],"_downGenes_WBID.txt ", partialPrefix, contrastNames[[grp]],"_down_tissue tissue\n"))
  cat(paste0("tea -d ../../../",pheDict," -q 0.05 -s ",partialPrefix,contrastNames[[grp]],"_downGenes_WBID.txt ", partialPrefix ,contrastNames[[grp]],"_down_phe phenotype\n"))
  cat(paste0("tea -d ../../../",goDict," -q 0.05 -s ",partialPrefix,contrastNames[[grp]],"_downGenes_WBID.txt ", partialPrefix, contrastNames[[grp]], "_down_go go\n"))
  sink()
}

sink(file=paste0(outPath,"/runTea.sh"),append=TRUE, type="output")
cat("cd ../../../\n")
sink()

system(paste0("grep -v '^>' ",outPath,"/runTea.sh > ",outPath,"/runTea1.sh"))
file.remove(paste0(outPath,"/runTea.sh"))
system(paste0("chmod +x ",outPath,"/runTea1.sh"))
system(paste0(outPath,"/runTea1.sh"),wait=F)




##########################-
##  Broadly Expressed genes Gernstein (2014)------
##########################-
broad<-read.csv(file=paste0(outPath,"/publicData/broadVregExpn_Gerstein2014.csv"),
        stringsAsFactors=T)

sigTables<-list()
bgCounts<-list()
for(grp in useContrasts){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],"_DESeq2_fullResults_p",padjVal,".rds"))

  sigTables[[paste0(grp)]]<-as.data.frame(
                        getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="both",
                        chr="all", nameChrCol="chr"))
  sigTables[[paste0(grp)]]$SMC<-grp
  bgCounts[[paste0(grp)]]<-as.data.frame(salmon)
  bgCounts[[paste0(grp)]]$SMC<-grp
}
sigGenes<-lapply(sigTables, "[", ,c("wormbaseID","sequenceID", "baseMean",
                                    "log2FoldChange","padj","SMC"))
bgCounts<-lapply(bgCounts,"[", ,c("wormbaseID","sequenceID", "baseMean",
                         "log2FoldChange","padj","SMC"))
lapply(sigGenes,dim)
lapply(bgCounts,dim)
#sigGenes<-lapply(sigGenes,na.omit)
sig<-do.call(rbind,sigGenes)
bgCounts<-do.call(rbind,bgCounts)

row.names(sig)<-NULL
sig<-inner_join(sig,broad,by="sequenceID")
sig$upVdown<-NA
sig$upVdown[sig$log2FoldChange>0]<-"up"
sig$upVdown[sig$log2FoldChange<0]<-"down"
sig$upVdown<-factor(sig$upVdown,levels=c("up","down"))

row.names(bgCounts)<-NULL
bgCounts<-inner_join(bgCounts,broad,by="sequenceID")
bgCounts$upVdown<-NA
bgCounts$upVdown[bgCounts$log2FoldChange>0]<-"up"
bgCounts$upVdown[bgCounts$log2FoldChange<0]<-"down"
bgCounts$upVdown<-factor(bgCounts$upVdown,levels=c("up","down"))


p1<-ggplot2::ggplot(sig, aes(x=category, y=abs(log2FoldChange), fill=upVdown))+
  geom_boxplot(notch=T,varwidth=F,outlier.shape=NA)+
  facet_wrap(~SMC)+ ylim(c(0,4))+theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("LFC of genes up/down regulated by domain type")


p2<-ggplot2::ggplot(sig, aes(x=category, y=log2(baseMean), fill=upVdown))+
  geom_boxplot(notch=T,varwidth=F,outlier.shape=NA)+
  facet_wrap(~SMC)+theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Mean expression of genes up/down regulated by domain type")


df<-sig %>% dplyr::group_by(SMC,category,upVdown) %>% dplyr::summarize(count=n())
dfbg<-bgCounts %>% dplyr::group_by(SMC,category) %>% dplyr::summarize(count=n())
df1<-left_join(dfbg,df,by=c("SMC","category"),suffix=c("_total",""))
df1$fraction<-df1$count/df1$count_total
p3<-ggplot2::ggplot(df1, aes(x=category, y=fraction, fill=upVdown))+
  geom_bar(stat="identity",position=position_dodge())+
  facet_wrap(~SMC)+theme_classic()+
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Fraction of genes") +
  ggtitle("Fraction of genes up/down regulated by domain type")

p<-ggpubr::ggarrange(p1,p2,p3,ncol=1,nrow=3)
ggplot2::ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"broadExpn_",
                                paste(useContrasts, collapse="_"),"_padj",
                                padjVal, "_lfc", lfcVal,".pdf"),
                plot=p, device="pdf",width=19,height=29,units="cm")


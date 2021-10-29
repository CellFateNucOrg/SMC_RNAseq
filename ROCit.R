library(ROCit)
library(ggplot2)
library(dplyr)
library(tidyr)


source("functions.R")
source("./variableSettings.R")

scriptName <- "ROCit"
print(scriptName)

if(filterData){
  fileNamePrefix<-filterPrefix
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
} else {
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
}
makeDirs(outPath,dirNameList=paste0(c("plots/"),
                                    paste0(dirname(fileNamePrefix),"/",
                                           scriptName)))



###########################
## compare samples
##########################


#######-
## venn diagrams------
#######-


## significantly changed genes
sigTables<-list()
for (grp in useContrasts){
  salmon<-readRDS(paste0(outPath,"/rds/",fileNamePrefix,contrastNames[[grp]],
                         "_DESeq2_fullResults_p",padjVal,".rds"))
  sigTables[[grp]]<-as.data.frame(getSignificantGenes(salmon, padj=padjVal, lfc=0,
                                                      namePadjCol="padj",
                                                      nameLfcCol="log2FoldChange",
                                                      direction="both",
                                                      chr="all", nameChrCol="chr"))
}
pdf(file=paste0(outPath,"/plots/",outputNamePrefix,"ROCit.pdf"), paper="a4r",
    height=8,width=11)
par(mfrow=c(2,3))
names(sigTables)
for (grp in useContrasts){
  tbl<-sigTables[[grp]]

  tbl$AvX<-ifelse(tbl$chr=="chrX","X","A")
  tbl$DC<-factor(ifelse(tbl$log2FoldChange>lfcVal & tbl$AvX=="X","DC","NDC"),
                 levels=c("NDC","DC"))

  measure <- measureit(score = abs(tbl$log2FoldChange), class = tbl$DC,
                       measure = c("ACC", "SENS", "FSCR"))
  names(measure)
  #> [1] "Cutoff" "Depth"  "TP"     "FP"     "TN"     "FN"     "ACC"    "SENS"
  #> [9] "FSCR"
  plot(measure$SENS~measure$Cutoff, type = "l",ylab="Measure")
  lines(measure$ACC~measure$Cutoff, col="blue")
  legend("bottomright",legend=c("SENS","ACC"),col=c("black","blue"),lty=c(1,1))
  crossPoint<-measure$Cutoff[which.min(abs(measure$ACC-measure$SENS))]
  maxPoint<-measure$Cutoff[which.max(measure$ACC+measure$SENS)]
  legend("topright",legend=c(paste0("Crossing point: ",round(crossPoint,2)),
                             paste0("Max SENS+ACC: ",round(maxPoint,2))))
  roc_empirical <- rocit(score = abs(tbl$log2FoldChange), class = tbl$DC,
                         negref = "NDC")
  summary(roc_empirical)
  names(roc_empirical)

  plot(roc_empirical, values = F)

  # KS plot-------

  rocitObj <- rocit(score = abs(tbl$log2FoldChange),
                 class = tbl$DC) #default: empirical
  kplot <- ksplot(rocitObj)
  legend("topright",legend=paste0("KS cutoff:",round(kplot$`KS Cutoff`,2)))
  print(kplot$`KS Cutoff`)
}
#print(kplot$`KS stat`)
dev.off()


library(readxl)
library(ggplot2)
library(EnhancedVolcano)
library(plyr)
library(dplyr)
library(ggpubr)
library(tidyr)


source("functions.R")
source("./variableSettings.R")

scriptName <- "compareGeneLists"
print(scriptName)

#chosenSubset=names(contrastNames) # for all complex contrasts
#chosenSubset=useContrasts[c(3,6,7,8)] #for most biologicall interesting contrasts
chosenSubset=useContrasts

#filterPrefix<-"p0.05_lfc0.5_filtChrAX/filtChrAX_"

if(filterData){
  fileNamePrefix<-filterPrefix
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
} else {
  outputNamePrefix<-gsub("\\/",paste0("/",scriptName,"/"),fileNamePrefix)
}

makeDirs(outPath,dirNameList=paste0(c("plots/"),
                                    paste0(dirname(fileNamePrefix),"/",
                                           scriptName)))


######################-
# Make combined table ----------------------------------------------
######################-
geneTable<-NULL
for (grp in chosenSubset){
  salmon<-as.data.frame(readRDS(paste0(outPath,"/rds/", fileNamePrefix, contrastNames[[grp]],
                                       "_DESeq2_fullResults_p",padjVal,".rds")))
  colnames(salmon)[colnames(salmon)=="baseMean"]<-paste0(grp,"_baseMean")
  colnames(salmon)[colnames(salmon)=="log2FoldChange"]<-paste0(grp,"_lfc")
  colnames(salmon)[colnames(salmon)=="padj"]<-paste0(grp,"_padj")
  if(is.null(geneTable)){
    geneTable<-as.data.frame(salmon)[,c("wormbaseID","chr", "start", "end",
                                        "strand", "publicID", "sequenceID",
                                      paste0(grp,"_baseMean"),
                                      paste0(grp,"_lfc"),
                                      paste0(grp,"_padj"))]
  } else {
    geneTable<-inner_join(geneTable,salmon[,c("wormbaseID","chr", "start",
                                              "end", "strand", "publicID",
                                              "sequenceID",
                                              paste0(grp,"_baseMean"),
                                              paste0(grp,"_lfc"),
                                              paste0(grp,"_padj"))],
                          by=c("wormbaseID","chr","start", "end",
                               "strand", "publicID", "sequenceID"))
  }
}



###################-
# cell cycle genes --------------------------------------------------------
###################-

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6404260/
#https://mcb.asm.org/content/24/6/2215
cellcycle<-read.delim(file=paste0(outPath,"/otherData/cellcycleGenes.txt"),header=T)

cellcycle<-join(cellcycle,geneTable,by="publicID")
cellcycle$publicID<-factor(cellcycle$publicID,levels=unique(cellcycle$publicID))
#cellcycle<-cellcycle[grep("^cy",cellcycle$publicID),]

lfcCols<-grep("_lfc",names(cellcycle))
df<-cbind(cellcycle[c("publicID","complex","role")],cellcycle[,lfcCols])
df <-df %>% tidyr::pivot_longer(cols=4:dim(df)[2],names_to="sample",values_to="Log2FoldChange")
df$sample<-gsub("_lfc","",df$sample)
df$sample<-factor(df$sample,levels=c("dpy26","kle2","scc1","coh1"))
colnames(df)[1]<-"gene"
#roleLab<-rle(df$role[df$sample==df$sample[1]])
#labPos<-cumsum(c(1,roleLab$lengths[-length(roleLab$lengths)]))#+(roleLab$lengths/2)
#labPos<-cumsum(rle(as.numeric(factor(df$role,levels=unique(df$role))))$lengths)

p1<-ggplot(df,aes(x=gene,y=Log2FoldChange)) +
  geom_bar(aes(fill=sample), stat="identity",position="dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  ggtitle("Cell cycle genes - log2 fold change") +
  geom_hline(yintercept=c(-0.5,0.5),linetype="dashed", color="grey60")#+
  #annotate("text", x=c(labPos), y = - 1, label =roleLab$values, size=2)

p1

padjCols<-grep("_padj",names(cellcycle))

df<-cbind(cellcycle$publicID,cellcycle[,padjCols])

df <-df %>% tidyr::gather(key="sample",value="padj",2:dim(df)[2])
df$sample<-gsub("_padj","",df$sample)
df$sample<-factor(df$sample,levels=c("dpy26","kle2","scc1","coh1"))
colnames(df)[1]<-"gene"

p2<-ggplot(df,aes(x=gene,y=-log(padj,base=10))) +
  geom_bar(aes(fill=sample), stat="identity",position="dodge") +
  theme_minimal() + ylab("-log10(P adjusted)") +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  ggtitle("Cell cycle genes - adjusted p value") +
  geom_hline(yintercept=-log(0.05,base=10),linetype="dashed", color="grey60")

p<-ggarrange(p1,p2,nrow=2)
p
ggsave(paste0(outPath,"/plots/",outputNamePrefix,"cellcycle.pdf"),height=19,width=29,units="cm",device="pdf")

#http://www.wormbook.org/chapters/www_cellcyclereguln/cellcyclereguln.html#sec2_1
#download.file("http://www.wormbook.org/chapters/www_cellcyclereguln/cellcyclefig1.jpg",
#              destfile=paste0(outPath,"/plots/cellcycle_cyclins_wormbook.jpg"))

###################-
# dna checkpoint genes --------------------------------------------------------
###################-

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1356337/

checkpoint<-read.delim(file=paste0(outPath,"/otherData/checkpointGenes.txt"),header=F)

colnames(checkpoint)[1]<-"publicID"

checkpoint<-join(checkpoint,geneTable,by="publicID")
checkpoint$publicID<-factor(checkpoint$publicID,levels=unique(checkpoint$publicID))


lfcCols<-grep("_lfc",names(checkpoint))
df<-cbind(checkpoint[c("publicID")],checkpoint[,lfcCols])
df <-df %>% gather(key="sample",value="Log2FoldChange",2:dim(df)[2])
df$sample<-gsub("_lfc","",df$sample)
colnames(df)[1]<-"gene"
#roleLab<-rle(df$role[df$sample==df$sample[1]])
#labPos<-cumsum(c(0,roleLab$lengths[-length(roleLab$lengths)]))+(roleLab$lengths/2)

p1<-ggplot(df,aes(x=gene,y=Log2FoldChange)) +
  geom_bar(aes(fill=sample), stat="identity",position="dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  ggtitle("DNA damage checkpoint - log2 fold change") +
  geom_hline(yintercept=c(-0.5,0.5),linetype="dashed", color="grey60") #+
#annotate("text", x=c(labPos), y = - 1, label =roleLab$values, size=2)


padjCols<-grep("_padj",names(checkpoint))

df<-cbind(checkpoint$publicID,checkpoint[,padjCols])

df <-df %>% gather(key="sample",value="padj",2:dim(df)[2])
df$sample<-gsub("_padj","",df$sample)
colnames(df)[1]<-"gene"

p2<-ggplot(df,aes(x=gene,y=-log(padj,base=10))) +
  geom_bar(aes(fill=sample), stat="identity",position="dodge") +
  theme_minimal() + ylab("-log10(P adjusted)") +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  ggtitle("DNA damage checkpoint - adjusted p value") +
  geom_hline(yintercept=-log(0.05,base=10),linetype="dashed", color="grey60")

p<-ggarrange(p1,p2,nrow=2)

ggsave(paste0(outPath,"/plots/",outputNamePrefix,"checkpoint.pdf"),height=19,width=19,units="cm",device="pdf")




###################-
# SMC complex genes --------------------------------------------------------
###################-

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1356337/

smc<-read.delim(file=paste0(outPath,"/otherData/SMCgenes.txt"),header=T)

#colnames(smc)[1]<-"publicID"

smc<-join(smc,geneTable,by="publicID")
smc$publicID<-factor(smc$publicID,levels=unique(smc$publicID))


lfcCols<-grep("_lfc",names(smc))
df<-cbind(smc[c("publicID")],smc[,lfcCols])
df <-df %>% gather(key="sample",value="Log2FoldChange",2:dim(df)[2])
df$sample<-gsub("_lfc","",df$sample)
colnames(df)[1]<-"gene"
#roleLab<-rle(df$role[df$sample==df$sample[1]])
#labPos<-cumsum(c(0,roleLab$lengths[-length(roleLab$lengths)]))+(roleLab$lengths/2)

p1<-ggplot(df,aes(x=gene,y=Log2FoldChange)) +
  geom_bar(aes(fill=sample), stat="identity",position="dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  ggtitle("SMC complex genes - log2 fold change") +
  geom_hline(yintercept=c(-0.5,0.5),linetype="dashed", color="grey60") #+
#annotate("text", x=c(labPos), y = - 1, label =roleLab$values, size=2)


padjCols<-grep("_padj",names(smc))

df<-cbind(smc$publicID,smc[,padjCols])

df <-df %>% gather(key="sample",value="padj",2:dim(df)[2])
df$sample<-gsub("_padj","",df$sample)
colnames(df)[1]<-"gene"

p2<-ggplot(df,aes(x=gene,y=-log(padj,base=10))) +
  geom_bar(aes(fill=sample), stat="identity",position="dodge") +
  theme_minimal() + ylab("-log10(P adjusted)") +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  ggtitle("SMC complex genes - adjusted p value") +
  geom_hline(yintercept=-log(0.05,base=10),linetype="dashed", color="grey60")

p<-ggarrange(p1,p2,nrow=2)

ggsave(paste0(outPath,"/plots/",outputNamePrefix,"smc.pdf"),height=19,width=29,units="cm",device="pdf")


#######################-
## Classic dosage compensation genes-------
#######################-

pubDC1<-readRDS(paste0(outPath,"/publicData/published_DCgr.rds"))
#add her-1
metadata<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))
her1<-metadata[grep("her-1",metadata$publicID),]
pubDC1<-c(her1,pubDC1)

if(filterData){
  # remove filtered genes
  idx<-pubDC1$wormbaseID %in% toFilter
  pubDC1<-pubDC1[!idx,]
}


pubDC<-join(data.frame(pubDC1),geneTable,by="wormbaseID")
pubDC$publicID<-ifelse(pubDC$publicID!="",pubDC$publicID,pubDC$sequenceID)
pubDC$publicID<-factor(pubDC$publicID,levels=unique(pubDC$publicID))


lfcCols<-grep("_lfc",names(pubDC))
df<-cbind(pubDC[c("publicID")],pubDC[,lfcCols])
df <-df %>% gather(key="sample",value="Log2FoldChange",2:dim(df)[2])
df$sample<-gsub("_lfc","",df$sample)
colnames(df)[1]<-"gene"
#roleLab<-rle(df$role[df$sample==df$sample[1]])
#labPos<-cumsum(c(0,roleLab$lengths[-length(roleLab$lengths)]))+(roleLab$lengths/2)

p1<-ggplot(df,aes(x=gene,y=Log2FoldChange)) +
  geom_bar(aes(fill=sample), stat="identity",position="dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  ggtitle("Dosage compensated genes - log2 fold change") +
  geom_hline(yintercept=c(-0.5,0.5),linetype="dashed", color="grey60") #+
#annotate("text", x=c(labPos), y = - 1, label =roleLab$values, size=2)


padjCols<-grep("_padj",names(pubDC))

df<-cbind(pubDC$publicID,pubDC[,padjCols])

df <-df %>% gather(key="sample",value="padj",2:dim(df)[2])
df$sample<-gsub("_padj","",df$sample)
colnames(df)[1]<-"gene"

p2<-ggplot(df,aes(x=gene,y=-log(padj,base=10))) +
  geom_bar(aes(fill=sample), stat="identity",position="dodge") +
  theme_minimal() + ylab("-log10(P adjusted)") +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  ggtitle("Dosage compensated genes - adjusted p value") +
  geom_hline(yintercept=-log(0.05,base=10),linetype="dashed", color="grey60")

p<-ggarrange(p1,p2,nrow=2)

ggsave(paste0(outPath,"/plots/",outputNamePrefix,"pubDC.pdf"),height=19,width=29,units="cm",device="pdf")



#######################-
## chromatin complexes-------
#######################-

chrom<-read_excel(paste0(outPath,"/publicData/chromatinComplexes.xlsx"))

metadata<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))

chrom<-left_join(chrom,data.frame(metadata),by="publicID")

if(filterData){
  # remove filtered genes
  idx<-chrom$wormbaseID %in% toFilter
  chrom<-chrom[!idx,]
}


chrom<-left_join(chrom,geneTable,by="wormbaseID")
chrom[,grep("\\.y$",colnames(chrom))]<-NULL
colnames(chrom)<-gsub("\\.x$","",colnames(chrom))

#chrom$publicID<-ifelse(chrom$publicID!="",chrom$publicID,chrom$sequenceID)
chrom$publicID<-factor(chrom$publicID,levels=chrom$publicID)


lfcCols<-grep("_lfc",names(chrom))
df<-cbind(chrom[c("publicID","complex")],chrom[,lfcCols])
df <-df %>% gather(key="sample",value="Log2FoldChange",3:dim(df)[2])
df$sample<-gsub("_lfc","",df$sample)
colnames(df)[1]<-"gene"
#roleLab<-rle(df$role[df$sample==df$sample[1]])
#labPos<-cumsum(c(0,roleLab$lengths[-length(roleLab$lengths)]))+(roleLab$lengths/2)

df$complex<-factor(df$complex,levels=unique(df$complex))
#dfann<-df %>% dplyr::group_by(complex) %>%
#  dplyr::summarise(count=n()/7) %>%
#  mutate(start=cumsum(c(1,count[-length(count)])),end=cumsum(count)) %>%
#  mutate(startGene=df$gene[start],endGene=df$gene[end])



p1<-ggplot(df,aes(x=gene,y=Log2FoldChange)) +
  geom_bar(aes(fill=sample), stat="identity",position="dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  ggtitle("Chromatin complexes - log2 fold change") +
  geom_hline(yintercept=c(-0.5,0.5),linetype="dashed", color="grey60") #+
#annotate("text", x=c(labPos), y = - 1, label =roleLab$values, size=2)

padjCols<-grep("_padj",names(chrom))

df<-cbind(chrom$publicID,chrom[,padjCols])

df <-df %>% gather(key="sample",value="padj",2:dim(df)[2])
df$sample<-gsub("_padj","",df$sample)
colnames(df)[1]<-"gene"

p2<-ggplot(df,aes(x=gene,y=-log(padj,base=10))) +
  geom_bar(aes(fill=sample), stat="identity",position="dodge") +
  theme_minimal() + ylab("-log10(P adjusted)") +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  ggtitle("Chromatin complexes - adjusted p value") +
  geom_hline(yintercept=-log(0.05,base=10),linetype="dashed", color="grey60")

p<-ggarrange(p1,p2,nrow=2)

ggsave(paste0(outPath,"/plots/",outputNamePrefix,"chromatinComplexes.pdf"),height=19,width=29,units="cm",device="pdf")


######################-
## relative abundance of subunits
######################-

# dds<-readRDS(paste0(outPath,"/rds/dds_object.rds"))
# smc<-read.delim(file=paste0(outPath,"/otherData/SMCgenes.txt"),header=T)
#
# metadata<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))
#
# smc1<-dplyr::left_join(smc,as.data.frame(metadata),by="publicID")
#
# smc1<-smc1[!is.na(smc1$wormbaseID),]
#
# #complex="Cohesin"
# smc1$baseMean<-NA
# idx<-match(smc1$wormbaseID,rowData(dds)$gene)
# smc1[!is.na(idx),"baseMean"]<-rowData(dds)$baseMean[na.omit(idx)]
#
# pdf(paste0(outPath,"/plots/baseMeanAbundanceSMC.pdf"))
# complexes=c("Cohesin","CondensinII","CondensinI","CondensinIDC","SDC")
# par(mfrow=c(3,2))
# for(complex in complexes){
#   hist(log2(rowData(dds)$baseMean),breaks=50,xlab="log2 baseMean",
#      main=paste(complex),col.main="#2222FFFF",col="grey90")
#   abline(v=median(log2(rowData(dds)$baseMean)),col="red")
#   rug(log2(smc1$baseMean[smc1[,complex]==1]),col="#2222FFFF",ticksize=1,lwd=1,
#       side=3)
# }
#
# hist(log2(rowData(dds)$baseMean),breaks=50,xlab="log2 baseMean",
#      main=paste("RNApol"),col.main="#2222FFFF")
# abline(v=median(log2(rowData(dds)$baseMean)),col="red")
# rug(log2(smc1$baseMean[smc1[,"SMC"]=="RNApol"]),col="#2222FFFF",ticksize=1,lwd=1)
#
# hist(log2(rowData(dds)$baseMean),breaks=50,xlab="log2 baseMean",
#      main=paste("topoisomerase"),col.main="#2222FFFF")
# rug(log2(smc1$baseMean[smc1[,"SMC"]=="topoisomerase"]),col="#2222FFFF",ticksize=1,lwd=1)
# dev.off()
#
# par(mfrow=c(1,1))

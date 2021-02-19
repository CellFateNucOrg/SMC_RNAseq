library(readxl)
library(ggVennDiagram)
library(ggplot2)
library(EnhancedVolcano)
library(dplyr)
library(plyr)

source("functions.R")

outPath="."
padjVal=0.05
lfcVal=0
plotPDFs=F
fileNamePrefix="noOsc_"

fileList<-read.table(paste0(outPath,"/fastqList.txt"), stringsAsFactors=F,
                     header=T)

# extract the strain variable
strain<-factor(as.character(unique(fileList$sampleName),
                            levels=c("366","382","775","784")))
SMC<-strain
levels(SMC)<-c("wt","dpy26cs","kle2cs","scc1cs")

controlGrp<-levels(SMC)[1] # control group
groupsOI<-levels(SMC)[-1]




######################-
# Make combined table ----------------------------------------------
######################-
geneTable<-NULL
for (grp in groupsOI){
  salmon<-as.data.frame(readRDS(paste0(outPath,"/rds/", fileNamePrefix, grp,
                                       "_DESeq2_fullResults.rds")))
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


lfcCols<-grep("_lfc",names(cellcycle))
df<-cbind(cellcycle[c("publicID")],cellcycle[,lfcCols])
df <-df %>% gather(key="sample",value="Log2FoldChange",2:dim(df)[2])
df$sample<-gsub("_lfc","",df$sample)
colnames(df)[1]<-"gene"
#roleLab<-rle(df$role[df$sample==df$sample[1]])
#labPos<-cumsum(c(0,roleLab$lengths[-length(roleLab$lengths)]))+(roleLab$lengths/2)

p1<-ggplot(df,aes(x=gene,y=Log2FoldChange)) +
  geom_bar(aes(fill=sample), stat="identity",position="dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  ggtitle("Cell cycle genes - log2 fold change") +
  geom_hline(yintercept=c(-0.5,0.5),linetype="dashed", color="grey60") #+
  #annotate("text", x=c(labPos), y = - 1, label =roleLab$values, size=2)


padjCols<-grep("_padj",names(cellcycle))

df<-cbind(cellcycle$publicID,cellcycle[,padjCols])

df <-df %>% gather(key="sample",value="padj",2:dim(df)[2])
df$sample<-gsub("_padj","",df$sample)
colnames(df)[1]<-"gene"

p2<-ggplot(df,aes(x=gene,y=-log(padj,base=10))) +
  geom_bar(aes(fill=sample), stat="identity",position="dodge") +
  theme_minimal() + ylab("-log10(P adjusted)") +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  ggtitle("Cell cycle genes - adjusted p value") +
  geom_hline(yintercept=-log(0.05,base=10),linetype="dashed", color="grey60")

p<-ggarrange(p1,p2,nrow=2)

ggsave(paste0(outPath,"/plots/",fileNamePrefix,"_cellcycle.pdf"),height=19,width=29,units="cm",device="pdf")



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

ggsave(paste0(outPath,"/plots/",fileNamePrefix,"_checkpoint.pdf"),height=19,width=19,units="cm",device="pdf")




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

ggsave(paste0(outPath,"/plots/",fileNamePrefix,"_smc.pdf"),height=19,width=29,units="cm",device="pdf")


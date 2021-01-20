outPath="."

if(! file.exists(paste0(outPath,"/wbGeneGRandRpts_WS275.rds"))) {
  source(paste0(outPath,"/createMetadataObj.R"))
}
metadata<-readRDS(paste0(outPath,"/wbGeneGRandRpts_WS275.rds"))
md<-readRDS(paste0(outPath,"metadataTbl_genes-rpts.rds"))


# aggregate htseq data
dataset="union_none"
files<-list.files(path=paste0(outPath,"/htseq"),pattern=paste0(dataset,".txt"),
                  full.names=T)

pdf(file=paste0(outPath,"/plots/BWA_",dataset,"_aggregatedVnorm.pdf"),
    paper="a4",height=19,width=11)
par(mfrow=c(4,2))
#f=files[1]
for (f in files){
  df<-read.delim(f,header=T)
  colnames(df)<-c("ID","count")
  rptRows<-grep("rpt",df$ID)
  rptdf<-df[rptRows,]
  idx<-match(rptdf$ID,metadata$ID)
  rptdf$rptfamID<-metadata$rptfamID[idx]

  rptFam<-rptdf %>% group_by(rptfamID) %>% summarise(count=sum(count))
  colnames(rptFam)<-c("ID","count")
  dffam<-rbind(df[-rptRows,],rptFam)
  outfile=gsub(".txt","_rptFam.txt",f)
  write.table(dffam,outfile,quote=F,col.names=F,row.names=F,sep="\t")

  #normalise by repeat family size
  rptFamNorm<-rptdf %>% group_by(rptfamID) %>% summarise(count=sum(count),
                                                     famSize=n())
  rptFamNorm$count<-ceiling(rptFamNorm$count/rptFamNorm$famSize)
  rptFamNorm$famSize<-NULL
  colnames(rptFamNorm)<-c("ID","count")
  dffamNorm<-rbind(df[-rptRows,],rptFamNorm)
  outfile=gsub(".txt","_rptFamNorm.txt",f)
  write.table(dffamNorm,outfile,quote=F,col.names=F,row.names=F,sep="\t")
  plot(rptFam$count~rptFamNorm$count,pch=16,col="#55555577",
       main=gsub(".txt","",basename(f)))
  abline(lm(rptFam$count~rptFamNorm$count),col="blue",lty=2)
}
dev.off()

# aggregate STAR data
dataset="ribo0"
files<-list.files(path=paste0(outPath,"/bamSTARrpts"),
                  pattern=paste0(dataset,"_ReadsPerGene.out.tab"), full.names=T)
f=files[1]
for (f in files){
  df<-read.delim(f,header=F)
  colnames(df)<-c("ID","Fwd","Rev","both")
  rptRows<-grep("rpt",df$ID)
  rptdf<-df[rptRows,]
  idx<-match(rptdf$ID,metadata$ID)
  rptdf$rptfamID<-metadata$rptfamID[idx]

  rptFam<-rptdf %>% group_by(rptfamID) %>% summarise(Fcount=sum(Fwd),
                                                     Rcount=sum(Rev),
                                                     bothCount=sum(both))
  colnames(rptFam)<-c("ID","Fwd","Rev","both")
  dffam<-rbind(df[-rptRows,],rptFam)
  outfile=gsub(".out.tab","_rptFam.tab",f)
  write.table(dffam,outfile,quote=F,col.names=F,row.names=F,sep="\t")
}


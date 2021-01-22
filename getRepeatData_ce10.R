library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce10)

genome<-Celegans
dfamVer="Dfam_2.0"

#args<-commandArgs(trailingOnly=TRUE)
#workDir<-args[1]
print(paste0("workDir is: ", workDir))
workDir="~/Documents/MeisterLab/otherPeopleProjects/Moushumi/SMC_RNAseq_repeats"
#workDir="."

dfamURL=paste0("https://www.dfam.org/releases/",dfamVer,"/ce10_dfam.nrph.hits.gz")
download.file(url=dfamURL, destfile=paste0(workDir,"/ce10_",dfamVer,".nrph.hits.gz"),method="wget")
system(paste0("gunzip ",workDir,"/ce10_",dfamVer,".nrph.hits.gz"))
system(paste0("grep -v '^#' ",workDir,"/ce10_",dfamVer,".nrph.hits > ",workDir,
              "/ce10_",dfamVer,".nrph.hits1"))
repeats_dfam <- read.delim(paste0(workDir,"/ce10_",dfamVer,".nrph.hits1"),
                           row.names=NULL, header=F)
repeats_dfam$V16<-NULL
colnames(repeats_dfam)<-unlist(strsplit(readLines(paste0(workDir,"/ce10_",dfamVer,".nrph.hits"),n=1),"\t"))

file.remove(paste0(workDir,"/ce10_",dfamVer,".nrph.hits1"))
file.remove(paste0(workDir,"/ce10_",dfamVer,".nrph.hits"))

repeats.gr <- GRanges(seqnames=repeats_dfam[,"#sequence name"],
                      ranges=IRanges(start=apply(repeats_dfam[,c("envelope start",
                                                          "envelope end")],1,FUN=min),
                                     end=apply(repeats_dfam[,c("envelope start",
                                                         "envelope end")],1,FUN=max)),
                      strand=repeats_dfam[,"strand of hit"])


repeats.gr$source<-dfamVer
repeats.gr$type<-"exon"

repeats.gr$ID<-NA
repeats.gr$Name<-repeats_dfam$`model name`
repeats.gr$Alias<-repeats_dfam$`model accession`
repeats.gr$Description<-repeats_dfam$`description of target sequence`
repeats.gr$ID<-paste0("rpt",1:length(repeats.gr),"_",
                           seqnames(repeats.gr),":",
                      start(repeats.gr),"-",end(repeats.gr))

saveRDS(repeats.gr,paste0(workDir,"/repeats_ce10_",dfamVer,"_nr.rds"))
export.gff3(repeats.gr, paste0(workDir,"/repeats_ce10_",dfamVer,"_nr.gff3"))


## bedgraph for Todor
repeats_bg<-repeats.gr
mcols(repeats_bg)<-NULL
repeats_bg$score<-1
repeats_bg$name<-repeats.gr$Name
#export.bedGraph(repeats_bg,paste0(workDir,"/repeats_ce10_dfam_nr.bedgraph"))
export.bed(repeats_bg,paste0(workDir,"/repeats_ce10_",dfamVer,"_nr.bed"))


# make wormbase formatted chr names for gff3
seqlevelsStyle(genome)<-"ensembl"
seqnames(genome)<-paste0("CHROMOSOME_",seqnames(genome))
seqlevelsStyle(repeats.gr)<-"ensembl"
seqlevels(repeats.gr)<-paste0("CHROMOSOME_",seqlevels(repeats.gr))
seqinfo(repeats.gr)<-seqinfo(genome)
genome(repeats.gr)<-NA
export.gff3(repeats.gr, paste0(workDir,"/repeats_WS220_",dfamVer,"_nr.gff3"))


# get repeat sequences
gff<-import(paste0(workDir,"/repeats_WS220_",dfamVer,"_nr.gff3"),format="gff3")

gtfdf<-data.frame(seqname=seqnames(gff),
                  source=gff$source,
                  feature=gff$type,
                  start=start(gff),
                  end=end(gff),
                  score=".",
                  strand=strand(gff),
                  frame=".",
                  attribute=paste0('gene_id "',gff$ID,'"; gene_name "',
                                   gff$Name,'";'))

write.table(gtfdf,paste0(workDir,"/repeats_WS220_",dfamVer,"_nr.gtf"),sep="\t",
            col.names=F,row.names=F, quote=F)


#seqlevelsStyle(genome)<-"ensembl"
rptSeq<-getSeq(genome,gff)

names(rptSeq)<-paste(gff$ID,gff$Name,gff$Alias,sep=";")

writeXStringSet(rptSeq, paste0(workDir,"/repeats_WS220_",dfamVer,"_nr.fa.gz"),
                append=FALSE, compress="gzip", format="fasta")


#### make chrom.sizes
chromSizes<-data.frame(seqnames=seqnames(genome),seqlengths=seqlengths(genome))
write.table(chromSizes,file=paste0(workDir,"/ws220.chrom.sizes"),col.names=F,
            row.names=F,sep="\t",quote=F)

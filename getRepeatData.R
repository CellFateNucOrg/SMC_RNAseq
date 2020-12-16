library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)

genome=Celegans
args<-commandArgs(trailingOnly=TRUE)
workDir<-args[1]
#workDir="~/Documents/MeisterLab/otherPeopleProjects/Moushumi/SMC_RNAseq_repeats"


dfamURL=paste0("https://www.dfam.org/releases/Dfam_3.3/annotations/ce10/ce10_dfam.nrph.hits.gz")
download.file(url=dfamURL, destfile=paste0(workDir,"/ce10_dfam.nrph.hits.gz"))
system(paste0("gunzip ",workDir,"/ce10_dfam.nrph.hits.gz"))
repeats_dfam <- read.delim(paste0(workDir,"/ce10_dfam.nrph.hits"))


repeats.gr <- GRanges(seqnames=as(repeats_dfam[,"X.seq_name"],"Rle"),
                      ranges=IRanges(start=apply(repeats_dfam[,c("ali.st",
                                                          "ali.en")],1,FUN=min),
                                     end=apply(repeats_dfam[,c("ali.st",
                                                         "ali.en")],1,FUN=max)),
                      strand=Rle(repeats_dfam[,"strand"]),
                      score="")

repeats.gr$source<-"Dfam_3.3"
repeats.gr$type<-"Transcript"

repeats.gr$ID<-paste0(seqnames(repeats.gr),":",
                      start(repeats.gr),"-",end(repeats.gr))
repeats.gr$Name<-repeats_dfam$family_name
repeats.gr$Alias<-repeats_dfam$family_acc

#export.gff3(repeats.gr, "repeats_ce10_dfam_nr.gff3")

liftoverURL="http://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz"
download.file(url=liftoverURL,destfile=paste0(workDir,"/ce10ToCe11.over.chain.gz"))
system(paste0("gunzip ",workDir,"/ce10ToCe11.over.chain.gz"))

chain <- import(paste0(workDir,"/ce10ToCe11.over.chain"))
repeats_ce11.gr <- unlist(liftOver(repeats.gr,chain))
seqlevelsStyle(repeats_ce11.gr)<-"ensembl"
export.gff3(repeats_ce11.gr, paste0(workDir,"/repeats_ce11_dfam_nr.gff3"))

#file.remove(paste0(workDir,"/ce10ToCe11.over.chain.gz"))
file.remove(paste0(workDir,"/ce10ToCe11.over.chain"))
#file.remove(paste0(workDir,"/ce10_dfam.nrph.hits.gz"))
file.remove(paste0(workDir,"/ce10_dfam.nrph.hits"))
#file.remove(paste0(workDir,"/repeats_ce10_dfam_nr.gff3"))


# get repeat sequences
gff<-import(paste0(workDir,"/repeats_ce11_dfam_nr.gff3"),format="gff3")

rptSeq<-getSeq(Celegans,gff)

names(rptSeq)<-paste(gff$ID,gff$Name,gff$Alias,sep=";")

writeXStringSet(rptSeq, paste0(workDir,"/repeats_ce11_dfam_nr.fa.gz"),
                append=FALSE, compress="gzip", format="fasta")


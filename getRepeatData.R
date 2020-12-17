library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)

genome<-Celegans

args<-commandArgs(trailingOnly=TRUE)
workDir<-args[1]
print(paste0("workDir is: ", workDir))
#workDir="~/Documents/MeisterLab/otherPeopleProjects/Moushumi/SMC_RNAseq_repeats"


dfamURL=paste0("https://www.dfam.org/releases/Dfam_3.3/annotations/ce10/ce10_dfam.nrph.hits.gz")
download.file(url=dfamURL, destfile=paste0(workDir,"/ce10_dfam.nrph.hits.gz"),method="wget")
system(paste0("gunzip ",workDir,"/ce10_dfam.nrph.hits.gz"))
repeats_dfam <- read.delim(paste0(workDir,"/ce10_dfam.nrph.hits"))


repeats.gr <- GRanges(seqnames=repeats_dfam[,"X.seq_name"],
                      ranges=IRanges(start=apply(repeats_dfam[,c("ali.st",
                                                          "ali.en")],1,FUN=min),
                                     end=apply(repeats_dfam[,c("ali.st",
                                                         "ali.en")],1,FUN=max)),
                      strand=repeats_dfam[,"strand"])

repeats.gr$source<-"Dfam_3.3"
repeats.gr$type<-"exon"

repeats.gr$ID<-NA
repeats.gr$Name<-repeats_dfam$family_name
repeats.gr$Alias<-repeats_dfam$family_acc

#export.gff3(repeats.gr, "repeats_ce10_dfam_nr.gff3")

liftoverURL="http://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz"
download.file(url=liftoverURL,destfile=paste0(workDir,"/ce10ToCe11.over.chain.gz"),method="wget")
system(paste0("gunzip ",workDir,"/ce10ToCe11.over.chain.gz"))

chain <- import(paste0(workDir,"/ce10ToCe11.over.chain"))
# liftover splits some granges, so need to fuse them again
repeats_ce11.gr <- unlist(range(liftOver(repeats.gr,chain)))
table(width(repeats_ce11.gr)-width(repeats.gr)) #checking the differences are minor
repeats_ce11.gr$ID<-paste0("rpt",1:length(repeats_ce11.gr),"_",
                           seqnames(repeats_ce11.gr),":",
                      start(repeats_ce11.gr),"-",end(repeats_ce11.gr))

## bedgraph for Todor
#repeats_bg<-repeats_ce11.gr
#mcols(repeats_bg)<-NULL
#repeats_bg$score<-1
#repeats_bg$name<-repeats_ce11.gr$Alias
#export.bedGraph(repeats_bg,paste0(workDir,"/repeats_ce11_dfam_nr.bedgraph"))
#export.bed(repeats_bg,paste0(workDir,"/repeats_ce11_dfam_nr.bed"))
#tileWidth=5000
#tiles<-unlist(tileGenome(seqinfo(Celegans),tilewidth=tileWidth))
#tiles$score<-countOverlaps(tiles,repeats_bg,ignore.strand=TRUE)
#export(tiles,paste0(workDir,"/repeats_ce11_dfam_nr_",tileWidth,".bw"),
#       format="bigwig")

# make wormbase formatted chr names for gff3
seqlevelsStyle(repeats_ce11.gr)<-"ensembl"
export.gff3(repeats_ce11.gr, paste0(workDir,"/repeats_ce11_dfam_nr.gff3"))

#file.remove(paste0(workDir,"/ce10ToCe11.over.chain.gz"))
file.remove(paste0(workDir,"/ce10ToCe11.over.chain"))
#file.remove(paste0(workDir,"/ce10_dfam.nrph.hits.gz"))
file.remove(paste0(workDir,"/ce10_dfam.nrph.hits"))
#file.remove(paste0(workDir,"/repeats_ce10_dfam_nr.gff3"))


# get repeat sequences
gff<-import(paste0(workDir,"/repeats_ce11_dfam_nr.gff3"),format="gff3")

seqlevelsStyle(genome)<-"ensembl"
rptSeq<-getSeq(genome,gff)

names(rptSeq)<-paste(gff$ID,gff$Name,gff$Alias,sep=";")

writeXStringSet(rptSeq, paste0(workDir,"/repeats_ce11_dfam_nr.fa.gz"),
                append=FALSE, compress="gzip", format="fasta")


#library(Organism.dplyr)
#library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
#library("TxDb.Celegans.UCSC.ce11.refGene")
#library("TxDb.Celegans.UCSC.ce11.ensGene")
library(rtracklayer)
library(dplyr)
#library(tximport)
#library(GenomicFeatures)

source("./variableSettings.R")

####
### some variables
#####

genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6], IRanges(start=1, end=seqlengths(Celegans)[1:6]))

wbseqinfo<-seqinfo(Celegans)
seqnames(wbseqinfo)<-c(gsub("chr","",seqnames(Celegans)))
seqnames(wbseqinfo)<-c(gsub("^M$","MtDNA",seqnames(wbseqinfo)))
genome(wbseqinfo)<-genomeVer
ce11seqinfo<-seqinfo(Celegans)

david<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/david_wbid2entrez_WS278.txt")


###############################################################-
### create metadata
###############################################################-

if(!file.exists(paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.",
                       genomeVer, ".annotations.sqlite"))){
  dir.create(paste0(genomeDir,"/annotations"),recursive=T)
  system(paste0("curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/gff/c_elegans.PRJNA13758.",
                genomeVer,".annotations.gff3.gz -o ",genomeDir,
                "/annotations/c_elegans.PRJNA13758.",genomeVer,
                ".annotations.gff3.gz"))

  system(paste0("gunzip ",genomeDir,"/annotations/c_elegans.PRJNA13758.",
                genomeVer,".annotations.gff3.gz"))
  si<-seqinfo(Celegans)
  genome(si)<-genomeVer
  seqnames(si)<-gsub("M","MtDNA",gsub("chr","",seqnames(si)))
  wstxdb<-makeTxDbFromGFF(file=paste0(genomeDir,
                                      "/annotations/c_elegans.PRJNA13758.",
                                      genomeVer,".annotations.gff3"),
                          format="gff3",organism="Caenorhabditis elegans",
                          chrominfo=si)

  saveDb(wstxdb,paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.", genomeVer,
                       ".annotations.sqlite"))
  file.remove(paste0(genomeDir,"/annotations/c_elegans.PRJNA13758.",
                     genomeVer, ".annotations.gff3"))
}

#gene ids
if(!file.exists(paste0(genomeDir,"/annotations/c_elegans.PRJNA13758.", genomeVer,
                       ".geneIDs.txt"))){
  system(paste0("curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/geneIDs/c_elegans.PRJNA13758.", genomeVer, ".geneIDs.txt.gz -o ", genomeDir,
                "/annotations/c_elegans.PRJNA13758.", genomeVer,
                ".geneIDs.txt.gz"))
  system(paste0("gunzip ",genomeDir,"/annotations/c_elegans.PRJNA13758.", genomeVer,
                ".geneIDs.txt.gz"))
}

# load a txdb of wormbase data and create a tx2gene object
txdb<-loadDb(paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.", genomeVer,
                    ".annotations.sqlite"))
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

columns(txdb)
keytypes(txdb)
TxptByGene<-transcriptsBy(txdb, by = "gene")
length(TxptByGene)

geneGR<-unlist(range(TxptByGene))
mcols(geneGR)$wormbase<-names(geneGR)
genedf<-as.data.frame(geneGR)


# download gene id data from simplemine: https://wormbase.org/tools/mine/simplemine.cgi
# for entrez ids, load wormbaseID column into https://david.ncifcrf.gov/conversion.jsp
#geneIDs<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/simplemine_WS278_geneID.txt")
geneIDs<-read.csv(paste0(genomeDir,"/annotations/c_elegans.PRJNA13758.", genomeVer,
       ".geneIDs.txt"),header=F)
names(geneIDs)<-c("organismID","wormbaseID","publicID","sequenceID","status","bioType")

metadata<-inner_join(geneIDs, genedf,by=c("wormbaseID"="wormbase")) %>%
  dplyr::select(wormbaseID,publicID,sequenceID,seqnames,start, end, strand) %>%
  collect %>% GenomicRanges::GRanges()

#names(mcols(metadata))<-c("wormbaseID","publicID","sequenceID")

i<-which(metadata$wormbaseID %in% david$From)
j<-match(metadata$wormbaseID[i],david$From)
metadata$entrezID<-NA
metadata$entrezID[i]<-david$To[j]

#seqinfo(metadata)<-wbseqinfo
seqlevelsStyle(metadata)<-"ucsc"
seqinfo(metadata)<-ce11seqinfo
metadata<-sort(metadata)

saveRDS(metadata,paste0(outPath,"/wbGeneGR_WS275.rds"))

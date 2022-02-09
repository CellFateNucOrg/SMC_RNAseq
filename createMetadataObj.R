library(AnnotationDbi)
library(BSgenome.Celegans.UCSC.ce11)
library(rtracklayer)
library(dplyr)
library(GenomicFeatures)

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



# import gtf file created with cufflinks from gff in indexGenomeTranscripts.sh
gtfFile<-paste0(genomeDir,"/annotations/c_elegans.PRJNA13758.",
                genomeVer,".annotations.gtf")

if(rnaType!="mRNA"){
  gtf<-import(gtfFile,format="gtf")
  seqlevelsStyle(gtf)<-"ucsc"
  seqinfo(gtf)<-seqinfo(Celegans)
  if(rnaType=="ncRNA" & !(file.exists(paste0(rnaType,"GR_WS275.rds")))){
    ncRNAurl<-paste0("ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.",genomeVer,".ncRNA_transcripts.fa.gz")
    if(!file.exists(gsub("\\.gz","",basename(ncRNAurl)))){
      download.file(ncRNAurl,destfile=basename(ncRNAurl))
      system(paste0("gunzip ", basename(ncRNAurl)))
    }
    system(paste0("echo sequenceID biotype wormbaseID publicID > ",rnaType,".txt"))
    system(paste0("grep '^>' ",gsub("\\.gz","",basename(ncRNAurl))," >> ",rnaType,".txt"))
    df<-read.delim(paste0(rnaType,".txt"),sep=" ", header=T, fill=T)
    df$sequenceID<-gsub("^>","",df$sequenceID)
    df$biotype<-gsub("^biotype=","",df$biotype)
    df$wormbaseID<-gsub("^gene=","",df$wormbaseID)
    df$publicID<-gsub("^locus=","",df$publicID)
    #table(df$biotype)
    # asRNA                    lincRNA                      miRNA
    # 104                        191                        458
    # miRNA_primary_transcript      nc_primary_transcript                        ncRNA
    # 10                        659                       7781
    # piRNA                  pre_miRNA                       rRNA
    # 15363                        260                         22
    # scRNA                     snoRNA                      snRNA
    # 1                        346                        129
    # transposable_element_ncRNA                       tRNA
    # 2                        634

    idx<-match(df$wormbaseID,gtf$gene_name)
    sum(is.na(idx))
    df<-df[!is.na(idx),] # remove 2 genes that are now considered transposons
    gr<-gtf[idx[!is.na(idx)]] # get relevant gr
    mcols(gr)<-df # put in clean metadata
    table(seqnames(gr))
#    chrI  chrII chrIII  chrIV   chrV   chrX   chrM
#    1293   1632   1125  16276   2312   3296     24
    gr<-sort(gr)
    saveRDS(gr,paste0(rnaType,"GR_WS275.rds"))
    file.remove(gsub("\\.gz","",basename(ncRNAurl)))
    file.remove(paste0(rnaType,".txt"))
  }


  if(rnaType=="pseudoRNA" & !(file.exists(paste0(rnaType,"GR_WS275.rds")))){
    ncRNAurl<-paste0("ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.",genomeVer,".pseudogenic_transcripts.fa.gz")
    if(!file.exists(gsub("\\.gz","",basename(ncRNAurl)))){
      download.file(ncRNAurl,destfile=basename(ncRNAurl))
      system(paste0("gunzip ", basename(ncRNAurl)))
    }
    system(paste0("echo sequenceID wormbaseID  biotype > ",rnaType,".txt"))
    system(paste0("grep '^>' ",gsub("\\.gz","",basename(ncRNAurl))," >> ",rnaType,".txt"))
    df<-read.delim(paste0(rnaType,".txt"),sep=" ", header=T, fill=T)
    df$sequenceID<-gsub("^>","",df$sequenceID)
    df$wormbaseID<-gsub("^gene=","",df$wormbaseID)
    df$biotype<-"pseudogene"
    #dim(df)
    #1902    3
    idx<-match(df$wormbaseID,gtf$gene_name)
    sum(is.na(idx))
    df<-df[!is.na(idx),]
    gr<-gtf[idx[!is.na(idx)]] # get relevant gr
    mcols(gr)<-df # put in clean metadata
    table(seqnames(gr))
    #chrI  chrII chrIII  chrIV   chrV   chrX   chrM
    #182    325     99    361    784    151      0
    gr<-sort(gr)
    saveRDS(gr,paste0(rnaType,"GR_WS275.rds"))
    file.remove(gsub("\\.gz","",basename(ncRNAurl)))
    file.remove(paste0(rnaType,".txt"))
  }

  if(rnaType=="tnRNA" & !(file.exists(paste0(rnaType,"GR_WS275.rds")))){
    ncRNAurl<-paste0("ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.",genomeVer,".transposon_transcripts.fa.gz")
    if(!file.exists(gsub("\\.gz","",basename(ncRNAurl)))){
      download.file(ncRNAurl,destfile=basename(ncRNAurl))
      system(paste0("gunzip ", basename(ncRNAurl)))
    }
    system(paste0("echo sequenceID biotype wormbaseID > ",rnaType,".txt"))
    system(paste0("grep '^>' ",gsub("\\.gz","",basename(ncRNAurl))," >> ",rnaType,".txt"))
    df<-read.delim(paste0(rnaType,".txt"),sep=" ", header=T, fill=T)
    df$sequenceID<-gsub("^>","",df$sequenceID)
    df$biotype<-gsub("^type=","",df$biotype)
    df$wormbaseID<-gsub("^gene=","",df$wormbaseID)
    #dim(df)
    #366    3
    idx<-match(df$wormbaseID,gtf$gene_name)
    sum(is.na(idx)) #
    df<-df[!is.na(idx),] # remove 2 genes that pseudo genes
    gr<-gtf[idx[!is.na(idx)]] # get relevant gr
    mcols(gr)<-df # put in clean metadata
    table(seqnames(gr))
    #chrI  chrII chrIII  chrIV   chrV   chrX   chrM
    #62     58     35     70     81     58      0
    gr<-sort(gr)
    saveRDS(gr,paste0(rnaType,"GR_WS275.rds"))
    file.remove(gsub("\\.gz","",basename(ncRNAurl)))
    file.remove(paste0(rnaType,".txt"))
  }
}

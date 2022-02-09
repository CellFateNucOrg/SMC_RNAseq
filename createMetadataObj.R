library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(dplyr)

## variables ###
source("./variableSettings.R")

genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)
david<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/david_wbid2entrez_WS278.txt")
dfam<-readRDS(paste0("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/SMC_RNAseq_repeats/repeats_",ucscVer,"_",dfamVer,"_nr.rds"))
ce11seqinfo<-seqinfo(Celegans)

# Create metadata object --------------------------------------------------
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
  # extract all gene lines
  system(paste0('grep "gene" ',genomeDir,'/annotations/c_elegans.PRJNA13758.', genomeVer,'.annotations.gff3 > ',genomeDir,'/annotations/geneLines_', genomeVer,'.gff3'))
  # extract all exon lines
  system(paste0('grep "exon" ',genomeDir,'/annotations/c_elegans.PRJNA13758.', genomeVer,'.annotations.gff3 > ',genomeDir,'/annotations/exonLines_', genomeVer,'.gff3'))
  file.remove(paste0(genomeDir,"/annotations/c_elegans.PRJNA13758.",
                     genomeVer, ".annotations.gff3"))
}


gff<-import(paste0(genomeDir,"/annotations/geneLines_",genomeVer,".gff3"))
gff<-gff[gff$type!="SAGE_tag" & gff$type!="exon",]
mcols(gff)<-mcols(gff)[,c("Name","sequence_name","locus","biotype")]
metadata<-gff[!is.na(gff$Name)]
names(mcols(metadata))<-c("wormbaseID","sequenceID","publicID","bioType")

# download gene id data from simplemine: https://wormbase.org/tools/mine/simplemine.cgi
# for entrez ids, load wormbaseID column into https://david.ncifcrf.gov/conversion.jsp
#geneIDs<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/simplemine_WS278_geneID.txt")
#david<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/david_wbid2entrez_WS278.txt")

i<-which(metadata$wormbaseID %in% david$From)
j<-match(metadata$wormbaseID[i],david$From)
metadata$entrezID<-NA
metadata$entrezID[i]<-david$To[j]

#seqinfo(metadata)<-wbseqinfo
seqlevelsStyle(metadata)<-"ucsc"
seqinfo(metadata)<-ce11seqinfo
metadata<-sort(metadata)

######## add repeat data
seqlevelsStyle(dfam)<-"ucsc"
seqinfo(dfam)<-seqinfo(metadata)

# give both objects the same column names
mcols(dfam)[,c("type","score","phase")]<-NULL
names(mcols(dfam))<-c("source","rptID","rptfamName","rptfamID","rptType")
mcols(dfam)[,names(mcols(metadata))]<-as.character(NA)
mcols(metadata)[,c("source","rptID","rptfamName","rptfamID","rptType")]<-as.character(NA)
metadata$source<-genomeVer
mcols(dfam)<-mcols(dfam)[,match(names(mcols(metadata)),names(mcols(dfam)))]
dfam$bioType<-"repeat"
# make sure column types are the same
mcols(dfam)<-sapply(mcols(dfam),as.character)
mcols(metadata)<-sapply(mcols(metadata),as.character)
#merge
md<-c(metadata,dfam)

# create universal ID column
md$ID<-md$wormbaseID
md$ID[is.na(md$wormbaseID)]<-md$rptID[is.na(md$wormbaseID)]

# mark repeats that overlap another gene and/or exon
exons<-import(paste0(genomeDir,'/annotations/exonLines_', genomeVer,'.gff3'),
              format="gff3")
exons<-exons[exons$type=="exon"]
mcols(exons)<-mcols(exons)[,c("source","type","Parent")]
seqlevelsStyle(exons)<-"ucsc"
rptRows<-which(md$bioType == "repeat")
toIgnore<-which(md$bioType %in% c("repeat","transposon","transposon_protein_coding",
                                  "transposon_pseudogene"))


olstart<-findOverlaps(md[rptRows,],md[-toIgnore,],minoverlap=10L,type="start")
olend<-findOverlaps(md[rptRows,],md[-toIgnore,],minoverlap=10L,type="end")
olin<-findOverlaps(md[rptRows,],md[-toIgnore,],minoverlap=10L,type="within")

exons<-exons[which(exons$source=="WormBase")]
olexons<-findOverlaps(md[rptRows,],exons,minoverlap=1L,type="any")
olexonsEqual<-findOverlaps(md[rptRows,],exons,minoverlap=1L,type="equal") # should be empty
length(unique(queryHits(olstart))) #39
length(unique(queryHits(olend))) #36
length(unique(queryHits(olin))) # 26113
length(unique(queryHits(olexons))) #6399

md$overlap<-NA
md$overlap[rptRows]<-"OL_none"
md$overlap[rptRows][unique(c(queryHits(olstart), queryHits(olend),
                             queryHits(olin)))]<-"OL_gene"
md$overlap[rptRows][unique(queryHits(olexons))]<-"OL_exon"

#OL_exon OL_gene OL_none
#6399   21551   55261       7.6% overlap exons and 25.9% overlap noncoding part of gene
saveRDS(md,paste0(outPath,"/wbGeneGRandRpts_",genomeVer,"_",dfamVer,".rds"))



## aggregated repeat families
dfamagg<-as.data.frame(dfam) %>% group_by(rptfamID) %>%
  summarise(famSize=n(),
            rptfamName=unique(rptfamName),
            rptfamID=unique(rptfamID),
            rptType=unique(rptType))

mddf<-as.data.frame(md)
mddf$famSize<-1

rptdf<-data.frame(matrix(ncol=ncol(mddf),nrow=nrow(dfamagg)))
colnames(rptdf)<-colnames(mddf)
rptdf[,colnames(dfamagg)]<-dfamagg
rptdf$ID<-rptdf$rptfamID
rptdf$bioType<-"repeatFamily"
mddf<-rbind(mddf,rptdf)

table(mddf$bioType)

saveRDS(mddf,paste0(outPath,"/metadataTbl_genes-rpts_",genomeVer,"_",dfamVer,".rds"))



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

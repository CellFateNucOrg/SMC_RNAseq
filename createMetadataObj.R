library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)

## variables ###
outPath="."
genomeVer="WS275"
dfamVer="Dfam_3.3"
ucscVer="ce11"

genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)
david<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/david_wbid2entrez_WS278.txt")
dfam<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/SMC_RNAseq_repeats/repeats_ce11_dfam_nr.rds")
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
dfam<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/SMC_RNAseq_repeats/repeats_ce11_dfam_nr.rds")

seqlevelsStyle(dfam)<-"ucsc"
seqinfo(dfam)<-seqinfo(metadata)

# give both objects the same column names
mcols(dfam)[,c("source","type","score","phase")]<-NULL
names(mcols(dfam))<-c("rptID","rptfamName","rptfamID","rptType")
mcols(dfam)[,names(mcols(metadata))]<-as.character(NA)
mcols(metadata)[,c("rptID","rptfamName","rptfamID","rptType")]<-as.character(NA)
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
length(unique(queryHits(olstart))) #18
length(unique(queryHits(olend))) #17
length(unique(queryHits(olin))) # 26278
length(unique(queryHits(olexons))) #4954

md$overlap<-NA
md$overlap[rptRows]<-"OL_none"
md$overlap[rptRows][unique(c(queryHits(olstart), queryHits(olend),
                             queryHits(olin)))]<-"OL_gene"
md$overlap[rptRows][unique(queryHits(olexons))]<-"OL_exon"


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


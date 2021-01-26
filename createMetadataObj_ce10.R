library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce10)

## variables ###
outPath="."
genomeVer="WS220"
dfamVer="Dfam_2.0"
ucscVer="ce10"

genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)
david<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/david_wbid2entrez_WS278.txt")
dfam<-readRDS(paste0("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/SMC_RNAseq_repeats/repeats_",ucscVer,"_",dfamVer,"_nr.rds"))
ce10seqinfo<-seqinfo(Celegans)

# Create metadata object --------------------------------------------------
###############################################################-
### create metadata
###############################################################-

if(!file.exists(paste0(genomeDir, "/annotations/c_elegans.",
                       genomeVer, ".annotations.sqlite"))){
  dir.create(paste0(genomeDir,"/annotations"),recursive=T)
  system(paste0("curl ftp://ftp.wormbase.org/pub/wormbase/releases/",genomeVer,
                "/species/c_elegans/c_elegans.",
                genomeVer,".annotations.gff3.gz -o ",genomeDir,
                "/annotations/c_elegans.",genomeVer,
                ".annotations.gff3.gz"))

  system(paste0("gunzip ",genomeDir,"/annotations/c_elegans.",
                genomeVer,".annotations.gff3.gz"))
  si<-seqinfo(Celegans)
  genome(si)<-genomeVer
  seqlevelsStyle(si)<-"ensembl"
  seqnames(si)<-paste0("CHROMOSOME_",seqnames(si))
  wstxdb<-makeTxDbFromGFF(file=paste0(genomeDir,
                                      "/annotations/c_elegans.",
                                      genomeVer,".annotations.gff3"),
                          format="gff3",organism="Caenorhabditis elegans",
                          chrominfo=si)

  saveDb(wstxdb,paste0(genomeDir, "/annotations/c_elegans.", genomeVer,
                       ".annotations.sqlite"))

  gff<-import(paste0(genomeDir,"/annotations/c_elegans.",
                genomeVer,".annotations.gff3"),format="gff3")
  # extract all gene lines
  system(paste0('grep "gene" ',genomeDir,'/annotations/c_elegans.', genomeVer,'.annotations.gff3 > ',genomeDir,'/annotations/geneLines_', genomeVer,'.gff3'))
  # extract all exon lines
  system(paste0('grep "exon" ',genomeDir,'/annotations/c_elegans.', genomeVer,'.annotations.gff3 > ',genomeDir,'/annotations/exonLines_', genomeVer,'.gff3'))
  file.remove(paste0(genomeDir,"/annotations/c_elegans.",
                     genomeVer, ".annotations.gff3"))

  system(paste0("curl ftp://ftp.wormbase.org/pub/wormbase/releases/",genomeVer,
                "/species/c_elegans/annotation/geneIDs.",genomeVer,".gz -o ",genomeDir,
                "/annotations/c_elegans.",genomeVer,
                ".geneIDs.txt.gz"))
  system(paste0("gunzip ",genomeDir,"/annotations/c_elegans.",genomeVer,
                ".geneIDs.txt.gz"))
}


gff<-import(paste0(genomeDir,"/annotations/geneLines_",genomeVer,".gff3"))
seqlevels(gff)<-seqlevels(gff)[c(1:4,6,7,5)]
si<-seqinfo(Celegans)
genome(si)<-genomeVer
seqlevelsStyle(si)<-"ensembl"
seqnames(si)<-paste0("CHROMOSOME_",seqnames(si))
seqinfo(gff)<-si

gff<-gff[gff$type!="SAGE_tag" & gff$type!="CDS",]
gff<-gff[gff$source!="landmark" & gff$source!="history",]
mcols(gff)<-mcols(gff)[,c("gene","ID","locus","source")]
names(mcols(gff))<-c("wormbaseID","sequenceID","publicID","bioType")
gff$sequenceID<-gsub("Gene:","",gff$sequenceID)
gff$sequenceID<-gsub("Pseudogene:","",gff$sequenceID)
ids<-data.frame(ID=gff$sequenceID)

# merge split pseudogenes
dups<-ids[duplicated(ids),]
gr<-gff[gff$sequenceID %in% dups,]
grl<-split(gr,gr$sequenceID)
gr<-unlist(range(grl))
mcols(gr)<-data.frame(wormbaseID=NA,sequenceID=names(gr),publicID=NA,
                      bioType="Pseudogene")
#gr$sequenceID<-rownames(gr)
gff<-c(gff[! gff$sequenceID %in% dups,],gr)

# find loci with two genes (same sequenceID but different WB and public IDs)
geneIDs<-read.csv(paste0(genomeDir,"/annotations/c_elegans.",genomeVer,
                         ".geneIDs.txt"), header=F,
                  na.strings=c("NA","NaN",""," "))
colnames(geneIDs)<-c("wormbaseID","publicID","sequenceID")
geneIDs<-geneIDs[!is.na(geneIDs$sequenceID),]
dups<-geneIDs[duplicated(geneIDs$sequenceID),]
geneIDs<-geneIDs[!(geneIDs$wormbaseID %in% dups$wormbaseID),]

# add gene ID data to gff
ids<-data.frame(ID=gff$sequenceID)
ids<-left_join(ids,geneIDs,by=c("ID"="sequenceID"),keep=T)
ids<-left_join(ids,geneIDs,by=c("ID"="wormbaseID"),keep=T)
xidx<-grep("WBGene",ids$wormbaseID.x)
names(gff)<-gff$sequenceID
mcols(gff)[xidx,c("wormbaseID","sequenceID","publicID")]<-ids[xidx,c("wormbaseID.x","sequenceID.x","publicID.x")]
yidx<-grep("WBGene",ids$wormbaseID.y)
mcols(gff)[yidx,c("wormbaseID","sequenceID","publicID")]<-ids[yidx,c("wormbaseID.y","sequenceID.y","publicID.y")]

#add both duplicates to gff
firstDup<-gff[gff$sequenceID %in% dups$sequenceID,]
secondDup<-firstDup
mcols(secondDup)[,c("wormbaseID","sequenceID","publicID")]<- dups[,c("wormbaseID","sequenceID","publicID")]

gff<-c(gff[!(gff$sequenceID %in% dups$sequenceID)],firstDup,secondDup)

## 20 remaining pseudogenes named after transcripts ad not genes:
unfused<-gff[is.na(gff$wormbaseID)]
names(unfused)<-gsub("[a|b]$","",names(unfused))
grl<-split(unfused,names(unfused))
gr<-unlist(range(grl))
idx<-match(names(gr),geneIDs$sequenceID)
mcols(gr)<-geneIDs[idx,c("wormbaseID","sequenceID","publicID")]
mcols(gr)$bioType<-"Pseudogene"

idx<-match(unfused$sequenceID,gff$sequenceID)
gff<-c(gff[-idx],gr)
names(gff)<-NULL
gff<-sort(gff)
saveRDS(gff,file=paste0(genomeDir,"/annotations/wbGeneAnnotationGR_",genomeVer,".rds"))


metadata<-gff

# download gene id data from simplemine: https://wormbase.org/tools/mine/simplemine.cgi
# for entrez ids, load wormbaseID column into https://david.ncifcrf.gov/conversion.jsp
#geneIDs<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/simplemine_WS278_geneID.txt")
#david<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/david_wbid2entrez_WS278.txt")

i<-which(metadata$wormbaseID %in% david$From)
j<-match(metadata$wormbaseID[i],david$From)
metadata$entrezID<-NA
metadata$entrezID[i]<-david$To[j]

#seqinfo(metadata)<-wbseqinfo
seqlevels(metadata)<-gsub("CHROMOSOME_","",seqlevels(metadata))
seqlevelsStyle(metadata)<-"ucsc"
seqinfo(metadata)<-ce10seqinfo
metadata<-sort(metadata)

######## add repeat data
dfam<-readRDS(paste0("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/SMC_RNAseq_repeats/repeats_",ucscVer,"_",dfamVer,"_nr.rds"))

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
mcols(exons)<-mcols(exons)[,c("source","type","Parent","ID","locus")]
seqlevels(exons)<-gsub("CHROMOSOME_","",seqlevels(exons))
seqlevelsStyle(exons)<-"ucsc"
seqlevels(exons)<-seqlevels(exons)[c(1:4,6,7,5)]
rptRows<-which(md$bioType == "repeat")
toIgnore<-which(md$bioType %in% c("repeat","rRNA"))


olstart<-findOverlaps(md[rptRows,],md[-toIgnore,],minoverlap=10L,type="start")
olend<-findOverlaps(md[rptRows,],md[-toIgnore,],minoverlap=10L,type="end")
olin<-findOverlaps(md[rptRows,],md[-toIgnore,],minoverlap=10L,type="within")

#exons<-exons[which(exons$source=="WormBase")]
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

saveRDS(mddf,paste0(outPath,"metadataTbl_genes-rpts_",genomeVer,"_",dfamVer,".rds"))


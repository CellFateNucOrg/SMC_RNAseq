## 2020-07-16
## Some published data sets were collected. A big list of dosage
## compensated and non dosage compenstated genes was manually copied
## from PDF of supplementary data (supl tables 4 & 5).
## A shorter list of "classical" dosage compenstaed and escaper genes was
## curated manually from figures or text in papers, assuming that if they
## are used as examples they are higher confidence (BUT they are all from
## embryos...):
## Wheeler 2016 just mentions some control genes in text. pcr in supl fig
## Csankovszki 2004 has a few genes in a figure first identified in Meyer 1986
## Kruesi 2016 in supl fig 1 to fig 7 took DC genes from Jans and tested
## 27 (?) of them for GROseq changes +- SDC-2, i took those with FC>1.5
## Jans 2009 has a few genes it confirmed by qPCR in supl table 6

library(magrittr)

outPath="."
genomeVer="WS275"
genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)

# txdb<-AnnotationDbi::loadDb(paste0(genomeDir,
#                                    "/annotations/c_elegans.PRJNA13758.",
#                                    genomeVer, ".annotations.sqlite"))
# #columns(txdb) # what kind of data is retrievable
# #keytypes(txdb)
# k <- keys(txdb, keytype = "GENEID")
# geneChr <- AnnotationDbi::select(txdb, k, columns=c("CDSCHROM"),
#                                  keytype="GENEID")


srcref <- Organism.dplyr::src_organism("TxDb.Celegans.UCSC.ce11.refGene")
metadata<-dplyr::inner_join(dplyr::tbl(srcref, "id"),
                                   dplyr::tbl(srcref, "ranges_gene")) %>%
  dplyr::select(wormbase, alias, genename, gene_chrom,
                gene_start, gene_end, gene_strand) %>%
  dplyr::collect() %>% GenomicRanges::GRanges()



#######################
## manually curated from papers
#######################

pubDC<-data.table::fread(input="published_DC.txt")
pubNDC<-data.table::fread(input="published_Xescapers.txt")


pubDCgr<-metadata[metadata$wormbase %in% pubDC$wbid]

mcols(pubDCgr)<-cbind(mcols(pubDCgr),pubDC[match(pubDCgr$wormbase,pubDC$wbid),])
pubDCgr$wbid<-NULL
pubDCgr$alias<-NULL
pubDCgr

saveRDS(pubDCgr,file="published_DCgr.rds")


pubNDCgr<-metadata[metadata$wormbase %in% pubNDC$wbid]

mcols(pubNDCgr)<-cbind(mcols(pubNDCgr),pubNDC[match(pubNDCgr$wormbase,pubNDC$wbid),])
pubNDCgr$wbid<-NULL
pubNDCgr$alias<-NULL
pubNDCgr

saveRDS(pubNDCgr,file="published_NDCgr.rds")


#######################
## Jans 2009
#######################

JansDC<-data.table::fread(input="Jans2009_DC_suplTable4.txt")
JansNDC<-data.table::fread(input="Jans2009_notDC_suplTable5.txt")


JansDCgr<-metadata[metadata$wormbase %in% JansDC$WormBaseId]

mcols(JansDCgr)<-cbind(mcols(JansDCgr),JansDC[match(JansDCgr$wormbase,JansDC$WormBaseId),])
JansDCgr$WormBaseId<-NULL
JansDCgr$alias<-NULL
JansDCgr

saveRDS(JansDCgr,file="Jans2009_DCgr.rds")


JansNDCgr<-metadata[metadata$wormbase %in% JansNDC$WormBaseId]

mcols(JansNDCgr)<-cbind(mcols(JansNDCgr),JansNDC[match(JansNDCgr$wormbase, JansNDC$WormBaseId),])
JansNDCgr$WormBaseId<-NULL
JansNDCgr$alias<-NULL
JansNDCgr

saveRDS(JansNDCgr,file="Jans2009_NDCgr.rds")




#######################
## Kramer 2015
#######################
kramerURL<-"https://doi.org/10.1371/journal.pgen.1005698.s011"
kramerFileName="Kramer_2015_PlotGen_S3file.xlsx"
#download.file(url=kramerURL,destfile=kramerFileName)

kramer<-readxl::read_excel(kramerFileName,col_types=c(rep("text",3),rep("numeric",30)))
names(kramer)
kramer<-kramer[,c(1:3,grep("_L3_",names(kramer)))]

kramergr<-metadata[metadata$wormbase %in% kramer$Gene_WB_ID]

mcols(kramergr)<-cbind(mcols(kramergr),kramer[match(kramergr$wormbase, kramer$Gene_WB_ID),])
kramergr$Gene_WB_ID<-NULL
kramergr$alias<-NULL

saveRDS(kramergr,file="kramer2015_gr.rds")
kramergr<-readRDS(file="kramer2015_gr.rds")


lfcVal<-0.5
padjVal<-0.05
idx<-!is.na(kramergr$dpy27_RNAi_L3_padj) &
  kramergr$dpy27_RNAi_L3_padj < padjVal &
  kramergr$dpy27_RNAi_L3_log2_fold_change > lfcVal &
  seqnames(kramergr)=="chrX"
kramerdpy27dc<-kramergr[idx,]
colIdx<-grep("(set)|(dpy21)|(mixed_sex)",colnames(mcols(kramerdpy27dc)))
mcols(kramerdpy27dc)[,colIdx]<-NULL
saveRDS(kramerdpy27dc,file=paste0("kramer2015_chrXup_dpy27_lfc",
                             formatC(lfcVal,format="e",digits=0),"_p",
                             formatC(padjVal,format="e",digits=0),
                             "_gr.rds"))


idx<-!is.na(kramergr$dpy21_mutant_L3_padj) &
  kramergr$dpy21_mutant_L3_padj < padjVal &
  kramergr$dpy21_mutant_L3_log2_fold_change > lfcVal &
  seqnames(kramergr)=="chrX"
kramerdpy21dc<-kramergr[idx,]
colIdx<-grep("(set)|(dpy27)",colnames(mcols(kramerdpy21dc)))
mcols(kramerdpy21dc)[,c(colIdx)]<-NULL
saveRDS(kramerdpy21dc,file=paste0("kramer2015_chrXup_dpy21_lfc",
                                  formatC(lfcVal,format="e",digits=0),"_p",
                                  formatC(padjVal,format="e",digits=0),
                                  "_gr.rds"))




###############-
#  Meeuse et al 2020 - Oscillating genes ----------------------------------
###############-
# Grosshans lab
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7370751/

meeuseURL<-"https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20209498&file=msb209498-sup-0003-DatasetEV1.xlsx"
meeuseFileName<-"MSB-16-e9498-s003.xlsx"

download.file(url=meeuseURL,destfile=meeuseFileName)

meeuse<-readxl::read_excel(meeuseFileName,col_types=c(rep("text",3),rep("numeric",2),"text"))
oscillating<-meeuse[meeuse$Class=="Osc",]

write.table(oscillating,file=paste0(outPath,"/oscillatingGenes.tsv"),row.names=F,
            col.names=T,quote=F,sep="\t")

file.remove(meeuseFileName)


###############-
#  Latorre et al 2015 - Oscillating genes ----------------------------------
###############-

latorreURL<-"http://genesdev.cshlp.org/content/suppl/2015/03/03/29.5.495.DC1/Supplemental_TableS7.xlsx"
latorreFileName<-"Supplemental_TableS7.xlsx"

download.file(url=latorreURL,destfile=latorreFileName)

latorre<-readxl::read_excel(latorreFileName,col_names=F)
colnames(latorre)<-"Osc_Latorre2015"

metadata<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))
sum(latorre$Osc_Latorre2015 %in% metadata$sequenceID)
#3235
length(latorre$Osc_Latorre2015)
#3269

idx<-match(latorre$Osc_Latorre2015,metadata$sequenceID)
latorre$wormbaseID<-metadata$wormbaseID[idx]
latorre<-latorre[!is.na(latorre$wormbaseID),]

write.table(latorre,file=paste0(outPath,"/oscillatingGenes_latorre.tsv"),row.names=F,
            col.names=T,quote=F,sep="\t")

file.remove(latorreFileName)

osc<-read.delim("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/SMC_RNAseq_filtCyc/oscillatingGenes.tsv")
sum(latorre$Osc_Latorre2015 %in% osc$SequenceName)
#2473
dim(latorre)
#3269
dim(osc)
#3739


###############-
# Garrigues 2019 - Heatshock L2 -------------------------------------------
###############-

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6927752/pdf/elife-51139.pdf
garriguesURL<-"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6927752/bin/elife-51139-supp2.xlsx"
garriguesFileName<-"elife-51139-supp2.xlsx"

download.file(url=garriguesURL,destfile=garriguesFileName)

garrigues<-readxl::read_excel(garriguesFileName,na="NA")
garrigues<-garrigues[! is.na(garrigues[,"P-adj"]),]

hsUP<-getSignificantGenes(garrigues, padj=0.05, lfc=1,
                          namePadjCol="P-adj",
                          nameLfcCol="log2(FC)", direction="gt")
hsUP
saveRDS(hsUP,file="hsUp_garrigues2019.rds")

hsDOWN<-getSignificantGenes(garrigues, padj=0.05, lfc=-1,
                          namePadjCol="P-adj",
                          nameLfcCol="log2(FC)", direction="lt")
hsDOWN
saveRDS(hsDOWN,file="hsDown_garrigues2019.rds")

file.remove(garriguesFileName)

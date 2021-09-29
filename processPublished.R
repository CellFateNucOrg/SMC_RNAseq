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
library(GenomicRanges)
library(eulerr)
library(lattice)
library(ggplot2)
library(DESeq2)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(readxl)

source("./variableSettings.R")
source("./functions.R")

if(!dir.exists(paste0(outPath,"/publicData"))) {
  dir.create(paste0(outPath,"/publicData"))
}

eulerLabelsType<-c("counts")
ce11seqinfo<-seqinfo(Celegans)

# txdb<-AnnotationDbi::loadDb(paste0(genomeDir,
#                                    "/annotations/c_elegans.PRJNA13758.",
#                                    genomeVer, ".annotations.sqlite"))
# #columns(txdb) # what kind of data is retrievable
# #keytypes(txdb)
# k <- keys(txdb, keytype = "GENEID")
# geneChr <- AnnotationDbi::select(txdb, k, columns=c("CDSCHROM"),
#                                  keytype="GENEID")


# srcref <- Organism.dplyr::src_organism("TxDb.Celegans.UCSC.ce11.refGene")
# metadata1<-dplyr::inner_join(dplyr::tbl(srcref, "id"),
#                                    dplyr::tbl(srcref, "ranges_gene")) %>%
#   dplyr::select(wormbase, alias, genename, gene_chrom,
#                 gene_start, gene_end, gene_strand) %>%
#   dplyr::collect() %>% GenomicRanges::GRanges()

metadata<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))


#######################-
## manually curated from papers------
#######################-


pubDC<-readxl::read_excel(paste0(outPath,"/publicData/DCgenes_published.xlsx"),sheet="DC")
pubNDC<-readxl::read_excel(paste0(outPath,"/publicData/DCgenes_published.xlsx"),sheet="Xescapers")


pubDCgr<-metadata[metadata$wormbaseID %in% pubDC$wbid]

mcols(pubDCgr)<-cbind(mcols(pubDCgr),pubDC[match(pubDCgr$wormbaseID,pubDC$wbid),])
pubDCgr$wbid<-NULL
pubDCgr$alias<-NULL
pubDCgr$geneid<-NULL
pubDCgr$publicid<-NULL
pubDCgr

saveRDS(pubDCgr,file=paste0(outPath,"/publicData/published_DCgr.rds"))


pubNDCgr<-metadata[metadata$wormbaseID %in% pubNDC$wbid]

mcols(pubNDCgr)<-cbind(mcols(pubNDCgr),pubNDC[match(pubNDCgr$wormbaseID,pubNDC$wbid),])
pubNDCgr$wbid<-NULL
pubNDCgr$alias<-NULL
pubNDCgr$geneid<-NULL
pubNDCgr$publicid<-NULL

saveRDS(pubNDCgr,file=paste0(outPath,"/publicData/published_NDCgr.rds"))


#######################-
## Jans 2009-----
#######################-

JansDC<-data.table::fread(input=paste0(outPath,"/publicData/Jans2009_DC_suplTable4.txt"))
JansNDC<-data.table::fread(input=paste0(outPath,"/publicData/Jans2009_notDC_suplTable5.txt"))


JansDCgr<-metadata[metadata$wormbaseID %in% JansDC$WormBaseId]

mcols(JansDCgr)<-cbind(mcols(JansDCgr),JansDC[match(JansDCgr$wormbaseID,JansDC$WormBaseId),])
JansDCgr$WormBaseId<-NULL
JansDCgr$alias<-NULL
JansDCgr$GeneName<-NULL
JansDCgr

saveRDS(JansDCgr,file=paste0(outPath,"/publicData/Jans2009_DCgr.rds"))


JansNDCgr<-metadata[metadata$wormbaseID %in% JansNDC$WormBaseId]

mcols(JansNDCgr)<-cbind(mcols(JansNDCgr),JansNDC[match(JansNDCgr$wormbaseID, JansNDC$WormBaseId),])
JansNDCgr$WormBaseId<-NULL
JansNDCgr$alias<-NULL
JansNDCgr$Gene_Name<-NULL
JansNDCgr

saveRDS(JansNDCgr,file=paste0(outPath,"/publicData/Jans2009_NDCgr.rds"))


#######################-
## Kramer 2015-----
#######################-
kramerURL<-"https://doi.org/10.1371/journal.pgen.1005698.s011"
kramerFileName="Kramer_2015_PlotGen_S3file.xlsx"

if(remakeFiles){
  file.remove(paste0(outPath,"/publicData/kramer2015_L3_gr.rds"))
}

if(! file.exists(paste0(outPath,"/publicData/kramer2015_L3_gr.rds"))){
  download.file(url=kramerURL,
                destfile=paste0(outPath,"/publicData/",kramerFileName))

  kramer<-readxl::read_excel(paste0(outPath,"/publicData/",kramerFileName),
                             col_types=c(rep("text",3),rep("numeric",30)))
  names(kramer)
  kramer<-kramer[,c(1:3,grep("_L3_",names(kramer)))]

  kramergr<-metadata[metadata$wormbaseID %in% kramer$Gene_WB_ID]

  mcols(kramergr)<-cbind(mcols(kramergr),kramer[match(kramergr$wormbaseID, kramer$Gene_WB_ID),])
  kramergr$Gene_WB_ID<-NULL
  kramergr$alias<-NULL
  kramergr$Sequence_Name_Gene<-NULL
  kramergr$Gene_Public_Name<-NULL

  saveRDS(kramergr,file=paste0(outPath,"/publicData/kramer2015_L3_gr.rds"))
  #file.remove(paste0(outPath,"/publicData/",kramerFileName))
}

kramergr<-readRDS(file=paste0(outPath,"/publicData/kramer2015_L3_gr.rds"))

localLFC<-0.5
localPadj<-0.05

idx<-!is.na(kramergr$dpy27_RNAi_L3_padj) &
  kramergr$dpy27_RNAi_L3_padj < localPadj &
  kramergr$dpy27_RNAi_L3_log2_fold_change > localLFC&
  seqnames(kramergr)=="chrX"
kramerdpy27dc<-kramergr[idx,]
colIdx<-grep("(set)|(dpy21)|(mixed_sex)",colnames(mcols(kramerdpy27dc)))
mcols(kramerdpy27dc)[,colIdx]<-NULL
saveRDS(kramerdpy27dc,file=paste0(outPath,"/publicData/kramer2015_chrXup_dpy27_lfc",
                                  localLFC,"_p", localPadj, "_gr.rds"))


idx<-!is.na(kramergr$dpy21_mutant_L3_padj) &
  kramergr$dpy21_mutant_L3_padj < localPadj &
  kramergr$dpy21_mutant_L3_log2_fold_change > localLFC &
  seqnames(kramergr)=="chrX"
kramerdpy21dc<-kramergr[idx,]
colIdx<-grep("(set)|(dpy27)",colnames(mcols(kramerdpy21dc)))
mcols(kramerdpy21dc)[,c(colIdx)]<-NULL
saveRDS(kramerdpy21dc,file=paste0(outPath,"/publicData/kramer2015_chrXup_dpy21_lfc",
                                  localLFC,"_p", localPadj, "_gr.rds"))


###############################-
## Jans 2009 vs Kramer-----
###############################-
kramer<-as.data.frame(readRDS(file=paste0(outPath,"/publicData/kramer2015_L3_gr.rds")))
names(kramer)
JansDC<-as.data.frame(readRDS(file=paste0(outPath,"/publicData/Jans2009_DCgr.rds")))
JansNDC<-as.data.frame(readRDS(file=paste0(outPath,"/publicData/Jans2009_NDCgr.rds")))

localPadj=0.05
localLFC=0.25
kramerdpy27dc<-getSignificantGenes(kramer, padj=localPadj, lfc=localLFC,
                                   namePadjCol="dpy27_RNAi_L3_padj",
                                   nameLfcCol="dpy27_RNAi_L3_log2_fold_change",
                                   direction="both",
                                   chr="all", nameChrCol="seqnames",
                                   outPath=outPath)
dim(kramerdpy27dc)
kramerdpy21dc<-getSignificantGenes(kramer, padj=localPadj, lfc=localLFC,
                                   namePadjCol="dpy21_mutant_L3_padj",
                                   nameLfcCol="dpy21_mutant_L3_log2_fold_change",
                                   direction="both",
                                   chr="all", nameChrCol="seqnames",
                                   outPath=outPath)

DC<-list(JansDC=JansDC$wormbaseID, dpy27=kramerdpy27dc$wormbaseID,
        dpy21=kramerdpy21dc$wormbaseID)
names(DC)<-c("JansDC", "dpy-27", "dpy-21")

pdf(file=paste0(outPath,
                "/publicData/venn_Jans2009vKramer2015_padj", localLFC,"_lfc",
                localPadj,".pdf"),width=5, height=10, paper="a4")
txtLabels<-list()
txtLabels[paste0("% ",names(DC)[2]," in ",names(DC)[1])]<-round(100*length(intersect(DC[[1]],DC[[2]]))/length(DC[[2]]),1)
txtLabels[paste0("% ",names(DC)[3]," in ",names(DC)[1])]<-round(100*length(intersect(DC[[1]],DC[[3]]))/length(DC[[3]]),1)

fit<-euler(DC)
percentages<-paste0(txtLabels[[1]], names(txtLabels)[1], " & ",
                    txtLabels[[2]], names(txtLabels)[2])
p1<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste0("Jans DC vs Kramer(2015): lfc>", localLFC,
                    ", padj<", localPadj,"\n",percentages), fontsize=8, y=0.7))
print(p1)




DC<-list(JansNDC=JansNDC$wormbaseID, dpy27=kramerdpy27dc$wormbaseID,
        dpy21=kramerdpy21dc$wormbaseID)
names(DC)<-c("JansNDC", "dpy-27", "dpy-21")
txtLabels<-list()
txtLabels[paste0("% ",names(DC)[2]," in ",names(DC)[1])]<-round(100*length(intersect(DC[[1]],DC[[2]]))/length(DC[[2]]),1)
txtLabels[paste0("% ",names(DC)[3]," in ",names(DC)[1])]<-round(100*length(intersect(DC[[1]],DC[[3]]))/length(DC[[3]]),1)

fit<-euler(DC)
percentages<-paste0(txtLabels[[1]], names(txtLabels)[1], " & ",
                    txtLabels[[2]], names(txtLabels)[2])
p2<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste0("Jans NDC vs Kramer(2015): lfc>", localLFC,
                  ", padj<",localPadj,"\n",percentages), fontsize=8, y=0.7))
print(p2)
dev.off()

# p<-ggpubr::ggarrange(p1,p2,ncol=2,nrow=1)
# ggplot2::ggsave(filename=paste0(outPath,
#                                 "/publicData/venn_Jans2009vKramer2015_padj",
#                                 localLFC,"_lfc", localPadj,".pdf"),
#                 plot=p, device="pdf",width=29,height=11,units="cm")




###############-
#  Meeuse et al 2020 - Oscillating genes -----------
###############-
# Developmental function and state transitions of a gene expression oscillator in Caenorhabditis elegan Meeuse...Grosshans Mol Syst Biol (2020)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7370751/

meeuseURL<-"https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20209498&file=msb209498-sup-0003-DatasetEV1.xlsx"
meeuseFileName<-"MSB-16-e9498-s003.xlsx"
if(remakeFiles){
  file.remove(paste0(outPath,"/publicData/oscillatingGenes.tsv"))
}
if(!file.exists(paste0(outPath,"/publicData/oscillatingGenes.tsv"))){
  download.file(url=meeuseURL,
                destfile=paste0(outPath,"/publicData/",meeuseFileName))

  meeuse<-readxl::read_excel(paste0(outPath,"/publicData/",meeuseFileName),
                             col_types=c(rep("text",3),rep("numeric",2),"text"))
  oscillating<-meeuse[meeuse$Class=="Osc",]
  names(oscillating)[1:3]<-c("wormbaseID","publicID","sequenceID")

  write.table(oscillating,file=paste0(outPath,"/publicData/oscillatingGenes.tsv"),
              row.names=F, col.names=T,quote=F,sep="\t")

  file.remove(paste0(outPath,"/publicData/",meeuseFileName))
}

###############-
#  Latorre et al 2015 - Oscillating genes ------------------
###############-
# The DREAM complex promotes gene body H2A.Z for target repression Latorre..Ahringer (2015)
# https://pubmed.ncbi.nlm.nih.gov/25737279/
latorreURL<-"http://genesdev.cshlp.org/content/suppl/2015/03/03/29.5.495.DC1/Supplemental_TableS7.xlsx"
latorreFileName<-"Supplemental_TableS7.xlsx"
if(remakeFiles){
  file.remove(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv"))
}

if(! file.exists(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv"))){
  download.file(url=latorreURL,destfile=paste0(outPath,"/publicData/",latorreFileName))

  latorre<-readxl::read_excel(paste0(outPath,"/publicData/",latorreFileName),
                              col_names=F)
  colnames(latorre)<-"Osc_Latorre2015"

  #metadata<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))
  sum(latorre$Osc_Latorre2015 %in% metadata$sequenceID)
  #3235
  length(latorre$Osc_Latorre2015)
  #3269

  idx<-match(latorre$Osc_Latorre2015,metadata$sequenceID)
  latorre$wormbaseID<-metadata$wormbaseID[idx]
  latorre<-latorre[!is.na(latorre$wormbaseID),]

  write.table(latorre,file=paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv"),
              row.names=F,col.names=T,quote=F,sep="\t")

  file.remove(paste0(outPath,"/publicData/",latorreFileName))
}

latorre<-read.delim(file=paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv"),
                     header=T,sep="\t")
osc<-read.delim(paste0(outPath,"/publicData/oscillatingGenes.tsv"))
sum(latorre$Osc_Latorre2015 %in% osc$SequenceName)
#2473
dim(latorre)
#3269
dim(osc)
#3739

OSC<-list(Meeuse=osc$wormbaseID,Latorre=latorre$wormbaseID)
length(unique(unlist(OSC))) #4522
fit<-euler(OSC)
pdf(paste0(outPath, "/publicData/venn_MeeuseVsLatorre.pdf"),width=5, height=10,
    paper="a4")
p1<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste0("Meeuse(2020) vs Latorre(2015) oscillating genes"),
                   fontsize=8, y=0.7))
print(p1)
dev.off()

# p<-ggVennDiagram(x) +
#   ggtitle(label=paste0("Meeuse(2020) vs Latorre(2015) oscillating genes"))
#
# ggplot2::ggsave(filename=paste0(outPath, "/publicData/venn_MeeuseVsLatorre.pdf"),
#                 plot=p, device="pdf",width=12,height=11,units="cm")




###############-
# Garrigues 2019 - Heatshock L2 -----------------
###############-

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6927752/pdf/elife-51139.pdf
garriguesURL<-"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6927752/bin/elife-51139-supp2.xlsx"
garriguesFileName<-"elife-51139-supp2.xlsx"
if(remakeFiles){
  file.remove(paste0(outPath,"/publicData/hsUp_garrigues2019.rds"))
}
if(! file.exists(paste0(outPath,"/publicData/hsUp_garrigues2019.rds"))){
  download.file(url=garriguesURL,destfile=paste0(outPath,"/publicData/",garriguesFileName))

  garrigues<-readxl::read_excel(paste0(outPath,"/publicData/",garriguesFileName),
                                na="NA")
  garrigues<-garrigues[! is.na(garrigues[,"P-adj"]),]

  idx<-match(c("WormBase.ID","gene.name"),names(garrigues))
  names(garrigues)[idx]<-c("wormbaseID","publicID")

  localPadj=0.05
  localLFC=1
  hsUP<-getSignificantGenes(garrigues, padj=localPadj, lfc=localLFC,
                            namePadjCol="P-adj",
                            nameLfcCol="log2(FC)", direction="gt")
  dim(hsUP)
  saveRDS(hsUP,file=paste0(outPath,"/publicData/hsUp_garrigues2019.rds"))

  hsDOWN<-getSignificantGenes(garrigues, padj=localPadj, lfc= -localLFC,
                              namePadjCol="P-adj",
                              nameLfcCol="log2(FC)", direction="lt")
  dim(hsDOWN)
  saveRDS(hsDOWN,file=paste0(outPath,"/publicData/hsDown_garrigues2019.rds"))

  file.remove(paste0(outPath,"/publicData/",garriguesFileName))
}




###############################-
## get germline-soma genes from Boeck-Waterston_GR2016-----
###############################-

#soma: JK1107(glp-1(q224)) mid-L4 30 h post-L1 stage larvae.
#Dissected gonads were from N2 (wild type) animals grown for 48 h at 20°C
#post-L1 stage larvae; approximately 200 gonads dissected and isolated from
#carcasses


tcFile="expressionTC_Boeck-Waterston_GR2016"
if(remakeFiles){
  file.remove(paste0(outPath,"/publicData/germlineSomaGenes_Boeck2016.csv"))
}

if(!file.exists(paste0(outPath,"/publicData/germlineSomaGenes_Boeck2016.csv"))) {
  link="https://genome.cshlp.org/content/suppl/2016/09/20/gr.202663.115.DC1/Supplemental_Table_S2.gz"
  download.file(link,paste0("publicData/",tcFile,".gz"))
  system(paste0("gunzip -S txt publicData/",tcFile,".gz"))

  tcData<-read.table(paste0("publicData/",tcFile),stringsAsFactors=F,header=T)
  glCols<-c("L4_counts","L4b_counts","L4JK1107soma_counts","L4JK1107soma.2_counts",
            "YA_counts","N2_Yad.1_counts","N2_Ad_gonad.1.RZLI_counts")

  #match(tcData$WormbaseName,metadata$sequenceID)
  tcData<-dplyr::left_join(tcData,as.data.frame(metadata),by=c("WormbaseName"="sequenceID"))
  tcData<-tcData[!is.na(tcData$wormbaseID),]
  tcData$WormbaseName<-NULL

  glCounts<-tcData[,glCols]
  row.names(glCounts)<-tcData$wormbaseID
  coldata<-data.frame(sampleNames=glCols,
                      stage=factor(c("L4","L4","L4","L4","YA","YA","YA")),
                      germline=factor(c("mix","mix","soma","soma","mix","mix",
                                        "gonad"),levels=c("mix","soma","gonad")))

  glCounts<-glCounts[!is.na(row.names(glCounts)),]


  dds <- DESeqDataSetFromMatrix(countData = glCounts,
                                colData = coldata,
                                design = ~ stage+germline)

  featureData <- data.frame(gene=rownames(glCounts))
  mcols(dds) <- DataFrame(mcols(dds), featureData)
  mcols(dds)

  #remove genes with few reads
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]

  # estimate parameters
  dds <- DESeq(dds)


  resLFCgermline<- lfcShrink(dds,contrast=c("germline","gonad","mix"),type="ashr")
  plotMA(resLFCgermline,main="gonad")
  summary(resLFCgermline)

  resLFCgermline<-resLFCgermline[!is.na(resLFCgermline$padj),]

  localLFC=0.5
  localPadj=0.05
  germline<-resLFCgermline[resLFCgermline$padj<localPadj & resLFCgermline$log2FoldChange>localLFC,]
  nongl<-resLFCgermline[resLFCgermline$padj<localPadj & resLFCgermline$log2FoldChange< -localLFC,]
  germline$wormbaseID<-rownames(germline)
  nongl$wormbaseID<-rownames(nongl)
  dim(germline) #2749
  dim(nongl) #4993


  resLFCsoma<- lfcShrink(dds,contrast=c("germline","soma","mix"),type="ashr")
  plotMA(resLFCsoma,main="soma")

  resLFCsoma<-resLFCsoma[!is.na(resLFCsoma$padj),]
  soma<-resLFCsoma[resLFCsoma$padj<localPadj & resLFCsoma$log2FoldChange>localLFC,]
  nonsoma<-resLFCsoma[resLFCsoma$padj<localPadj & resLFCsoma$log2FoldChange< -localLFC,]
  soma$wormbaseID<-rownames(soma)
  nonsoma$wormbaseID<-rownames(nonsoma)
  dim(soma) #2127
  dim(nonsoma) #3144
  sum(soma$wormbaseID %in% germline$wormbaseID) #49
  sum(nonsoma$wormbaseID %in% germline$wormbaseID) #939
  sum(germline$wormbaseID %in% soma$wormbaseID) #49
  sum(nongl$wormbaseID %in% soma$wormbaseID) #1293
  sum(nongl$wormbaseID %in% nonsoma$wormbaseID) #449

  #germline<-germline[!(germline$wormbaseID %in% soma$wormbaseID),]
  #soma<-soma[!(soma$wormbaseID %in% germline$wormbaseID),]

  glvSoma<-data.frame(wormbaseID=c(germline$wormbaseID,nongl$wormbaseID,soma$wormbaseID,
                                   nonsoma$wormbaseID),
                      germline=c(rep("gonadYA",nrow(germline)),
                                 rep("somaYA",nrow(nongl)),
                                 rep("somaL4", nrow(soma)),
                                 rep("germlineL4",nrow(nonsoma))))
  write.csv(glvSoma,paste0(outPath,"/publicData/germlineSomaGenes_Boeck2016.csv"),
            row.names=F)

  file.remove(paste0("publicData/",tcFile))
}


###############################-
## get germline-soma genes from Reinke_DEV2004-----
###############################-
# https://journals.biologists.com/dev/article/131/2/311/42475/Genome-wide-germline-enriched-and-sex-biased
# Genome-wide germline-enriched and sex-biased expression profiles in Caenorhabditis elegans
# Valerie Reinke ,  Inigo San Gil ,  Samuel Ward ,  Keith Kazmer
# Development (2004) 131 (2): 311–323.

glFile="germlineVsoma_Reinke_Dev2004"

if(remakeFiles){
  file.remove(paste0(outPath,"/publicData/germlineGenes_Reinke2004.csv"))
}

if (!file.exists(paste0(outPath,"/publicData/germlineGenes_Reinke2004.csv"))) {
  #link2="http://dev.biologists.org/highwire/filestream/1201187/field_highwire_adjunct_files/0/Data_S1.zip"
  link2="https://cob.silverchair-cdn.com/cob/content_public/journal/dev/131/2/10.1242_dev.00914/5/dev00914-sup-data_s1.zip?Expires=1632688021&Signature=XWm42IDBQhktZa-HFkbXCUf9tsRusTi~1T0MCyfi3smx69Yl7BcM078~4ObJOKugz6Ojxb5aNlwT94Mf09p-5A7FvUnT0vCAf8xJGZ3Byst1eIIBR4v4M6Smb487D0dh2DbZCsWmogyh6XZNTOss0bOrFcxJ~NUg9fmVR59yYLBUiHuMiEbfHyVd3ACy7bRlcwiXRSsyI4NaWDHaPluFgqgMlFEsibvFhT3aKG3BWs2IbUTearqfls5o2mwxP0GFY5ZAaR973-K7qtNY7U2vKH-cdoO2Aec~OJOKh42VlUaE1F5T8pJ5H2mkT6ZXOK3BPTO7Gmz33vmGzta4nduogg__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA"
  download.file(link2,paste0("publicData/",glFile,".zip"))
  system(paste0("rm -rf ",outPath,"/publicData/",glFile))
  system(paste0("unzip ",outPath,"/publicData/",glFile,".zip -d ",outPath,"/publicData/",glFile))

  glData<-read.delim(paste0("publicData/",glFile,"/Fig1\ I\ wt\ vs\ glp4\ enriched\ genes.txt"),stringsAsFactors=F,header=T,sep="\t")
  glData<-glData[,c("WormbaseID","exclusive.category")]

  table(glData$exclusive.category)

  glData<-glData[glData$exclusive.category %in% c("herm intrinsic", "herm oocyte",
                                                  "herm sperm","shared intrinsic",
                                                  "shared oocyte","shared sperm"),]

  glData<-dplyr::left_join(glData, as.data.frame(metadata),by=c("WormbaseID"="sequenceID"))
  names(glData)[1]<-"sequenceID"
  glData<-glData[!is.na(glData$wormbaseID),]

  write.csv(glData,paste0(outPath,"/publicData/germlineGenes_Reinke2004.csv"),
            row.names=F)

  file.remove(paste0(outPath,"/publicData/",glFile,"/Fig1\ I\ wt\ vs\ glp4\ enriched\ genes.txt"))
  file.remove(paste0(outPath,"/publicData/",glFile,".zip"))
  system(paste0("rm -rf ",outPath,"/publicData/",glFile))
}



######################-
## compare germline datasets-------
######################-
pdf(paste0(outPath, "/publicData/venn_ReinkeVsBoeck.pdf"),width=5, height=10,
    paper="a4")

glData<-read.csv(paste0(outPath,"/publicData/germlineGenes_Reinke2004.csv"))
glvSoma<-read.csv(paste0(outPath,"/publicData/germlineSomaGenes_Boeck2016.csv"))
germsoma<-list(germline=glData$wormbaseID,
        germlineL4=glvSoma$wormbaseID[glvSoma$germline=="germlineL4"],
        somaL4=glvSoma$wormbaseID[glvSoma$germline=="somaL4"])

fit<-euler(germsoma)
totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
                sum(fit$original.values[grep(x, names(fit$original.values))]))}),
                 collapse="  ")
p1<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste0("Reinke(2004) vs Boeck(2016) L4 soma\n",totalSums),
                   fontsize=8))
print(p1)
#p1<-ggVennDiagram(germsoma) + ggtitle(label=paste0("Reinke(2004) vs Boeck(2016) L4 soma"))

germsoma<-list(germline=glData$wormbaseID,
        gonadYA=glvSoma$wormbaseID[glvSoma$germline=="gonadYA"],
        somaYA=glvSoma$wormbaseID[glvSoma$germline=="somaYA"])

fit<-euler(germsoma)
totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
        sum(fit$original.values[grep(x, names(fit$original.values))]))}),
                 collapse="  ")
p2<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste0("Reinke(2004) vs Boeck(2016) YA gonad\n",totalSums),
                   fontsize=8))
print(p2)
#p2<-ggVennDiagram(germsoma) + ggtitle(label=paste0("Reinke(2004) vs Boeck(2016) YA gonad"))

germsoma<-list(germlineL4=glvSoma$wormbaseID[glvSoma$germline=="germlineL4"],
        gonadYA=glvSoma$wormbaseID[glvSoma$germline=="gonadYA"],
        somaL4=glvSoma$wormbaseID[glvSoma$germline=="somaL4"],
        somaYA=glvSoma$wormbaseID[glvSoma$germline=="somaYA"])
fit<-euler(germsoma)
totalSums<-paste(lapply(row.names(fit$ellipses), function(x){paste(x,
             sum(fit$original.values[grep(x, names(fit$original.values))]))}),
                 collapse="  ")
p3<-plot(fit, quantities=list(type=eulerLabelsType),
         main=list(label=paste0("Boeck(2016) L4 vs YA soma/germline\n",totalSums),
                   fontsize=8))
print(p3)
#p3<-ggVennDiagram(germsoma) + ggtitle(label=paste0("Boeck(2016) L4 vs YA soma/germline"))
dev.off()
#p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
#ggplot2::ggsave(filename=paste0(outPath, "/publicData/venn_ReinkeVsBoeck.pdf"),
#                plot=p, device="pdf",width=29,height=11,units="cm")


# aging genes:
# LOF of daf-2 (insulin receptor) increases lifespan
# LOF of age-1 (PI3K) increases lifespan
# LOF of daf-16 (TF repressed by age-1) decreases lifespan
# LOF of elt-3 causes a decrease in lifespan
# elt-1 activates elt-3. elt-5 and elt-6 repress elt-3
# LOF of elt-5&elt-6 cause increase in lifespan

#######################-
## Aging microarrays Budovskaya (2008) (Stuart Kim lab)------
#######################-
# https://www.cell.com/fulltext/S0092-8674%2808%2900707-1

# Probably same data as:  https://www.sciencedirect.com/science/article/pii/S0960982202011466?via%3Dihub
#Supplemental Table 1. Gene expression levels during aging in C. elegans.
#http://cmgm.stanford.edu/~kimlab/aging/suptable1.txt

#also see:
#https://science.sciencemag.org/content/293/5537/2087.full?ijkey=MsA0e.Sfl1Wpw&keytype=ref&siteid=sci
xyTopoURL<-"http://cmgm.stanford.edu/~kimlab/topomap/worm3.txt"


# Document S4. Table S3: Aging Microarray Data
agingMAurl<-"https://www.cell.com/cms/10.1016/j.cell.2008.05.044/attachment/424c8398-9071-4f6f-a500-be367eb342ca/mmc4.xls"

# Document S5. Table S4: 1254 Age-Regulated Genes
ageRegulatedURL<-"https://www.cell.com/cms/10.1016/j.cell.2008.05.044/attachment/bf38a101-cfb0-48cc-b421-b108a59bd267/mmc5.xls"

if(remakeFiles){
  file.remove(paste0(outPath,"/publicData/AgeRegulated_Budovskaya2008.csv"))
}

if(!file.exists(paste0(outPath,"/publicData/AgeRegulated_Budovskaya2008.csv"))){
  download.file(ageRegulatedURL,paste0(outPath,"/publicData/mmc5_Budovskaya.xls"))
  ageReg<-readxl::read_excel(paste0(outPath,"/publicData/mmc5_Budovskaya.xls"))
  colnames(ageReg)<-c("sequenceID") #1244 but only 1011 overlap with ws275 gene names
  idx<-ageReg$sequenceID %in% metadata$sequenceID
  ageReg<-ageReg[idx,]
  ageReg<-left_join(ageReg,as.data.frame(metadata),by=c("sequenceID"))

  download.file(agingMAurl,paste0(outPath,"/publicData/mmc4_Budovskaya.xls"))
  tcMA<-readxl::read_excel(paste0(outPath,"/publicData/mmc4_Budovskaya.xls"),skip=4)
  colnames(tcMA)[1]<-"sequenceID"
  idx<-tcMA$sequenceID %in% metadata$sequenceID
  tcMA<-tcMA[idx,]
  tcMA<-left_join(tcMA,as.data.frame(metadata),by=c("sequenceID"))
  write.csv(tcMA,paste0(outPath,"/publicData/AgingTC_Budovskaya2008.csv"),
            row.names=F,quote=F)

  ageReg<-left_join(ageReg,tcMA,by=c("sequenceID"="gene"))
  write.csv(ageReg,paste0(outPath,"/publicData/AgeRegulated_Budovskaya2008.csv"),
            row.names=F,quote=F)
  file.remove(paste0(outPath,"/publicData/mmc5_Budovskaya.xls"))
  file.remove(paste0(outPath,"/publicData/mmc4_Budovskaya.xls"))
}

# Document S6. Table S5: age-1 Microarray Data
age1MAurl<-"https://www.cell.com/cms/10.1016/j.cell.2008.05.044/attachment/b1a6e947-10b7-438b-8a3a-690c71964854/mmc6.xls"
if(!file.exists(paste0(outPath,"/publicData/age1MA_Budovskaya2008.csv"))){
  download.file(age1MAurl,paste0(outPath,"/publicData/mmc6_Budovskaya.xls"))
  # there is some in the file. need to open manually and save as xlsx
  age1MA<-readxl::read_excel(paste0(outPath,"/publicData/mmc6_Budovskaya.xlsx"),skip=9)
  # keep only final summary columns
  age1MA<-age1MA[,c(1,2,19,20,21,22,23)]
  # rename columns
  colnames(age1MA)<-c("sequenceID","description","ratio","std","df","tval","pval")
  #idx<-age1MA$sequenceID %in% metadata$sequenceID #18556 genes down to 14645
  #age1MA<-age1MA[idx,]
  age1MA<-inner_join(age1MA,as.data.frame(metadata),by="sequenceID")
  write.table(age1MA,paste0(outPath,"/publicData/age1MA_Budovskaya2008.tsv"),
            row.names=F,quote=F,sep="\t")
  file.remove(paste0(outPath,"/publicData/mmc6_Budovskaya.xls"))
}


# Document S7. Table S6: daf-16(m26) Microarray Data
daf16MAurl<-"https://www.cell.com/cms/10.1016/j.cell.2008.05.044/attachment/e172c477-a9c6-4c52-bb8f-fa5211ad4cee/mmc7.xls"
if(!file.exists(paste0(outPath,"/publicData/daf16MA_Budovskaya2008.csv"))){
  download.file(daf16MAurl,paste0(outPath,"/publicData/mmc7_Budovskaya.xls"))
  daf16MA<-readxl::read_excel(paste0(outPath,"/publicData/mmc7_Budovskaya.xls"),skip=9)
  # keep only final summary columns
  daf16MA<-daf16MA[,c(1,2,19,20,21,22,23)]
  # rename columns
  colnames(daf16MA)<-c("sequenceID","description","ratio","std","df","tval","pval")
  #idx<-daf16MA$sequenceID %in% metadata$sequenceID #18556 genes down to 14645
  #daf16MA<-daf16MA[idx,]
  daf16MA<-inner_join(daf16MA,as.data.frame(metadata),by="sequenceID")
  write.table(daf16MA,paste0(outPath,"/publicData/daf16MA_Budovskaya2008.tsv"),
            row.names=F,quote=F,sep="\t")
  file.remove(paste0(outPath,"/publicData/mmc7_Budovskaya.xls"))
}



#Lund (2002) paper:
#  https://pubmed.ncbi.nlm.nih.gov/12372248/
# aging timecourse data hard to extract



#######################-
## Aging microarrays Murphy(2003)------
#######################-

#Murphy (2003) paper
#https://www.nature.com/articles/nature01789#MOESM1
# daf-2,  daf-2/daf-16 MA
# I want classI genes upregulated by daf-2 (RNAi) and downregulated by daf-16 RNAi (=pro longevity)
murphyURL<-"https://static-content.springer.com/esm/art%3A10.1038%2Fnature01789/MediaObjects/41586_2003_BFnature01789_MOESM2_ESM.xls"

if(!file.exists(paste0(outPath,"/publicData/agingClassI_Murphy2003.csv"))){
  download.file(murphyURL,paste0(outPath,"/publicData/41586_2003_BFnature01789_MOESM2_ESM.xls"))
  # there is some in the file. need to open manually and save as xlsx
  ageClassI<-readxl::read_excel(paste0(outPath,"/publicData/41586_2003_BFnature01789_MOESM2_ESM.xls"),
                                sheet=1,skip=1,col_names=c("sequenceID","description"))
  ageClassII<-readxl::read_excel(paste0(outPath,"/publicData/41586_2003_BFnature01789_MOESM2_ESM.xls"),
                                sheet=2,skip=1,col_names=c("sequenceID","description"))
  ageClassI$sequenceID<-gsub("\\*+","",ageClassI$sequenceID)
  ageClassII$sequenceID<-gsub("\\*+","",ageClassII$sequenceID)
  ageClassI<-inner_join(ageClassI,as.data.frame(metadata),by="sequenceID")
  ageClassII<-inner_join(ageClassII,as.data.frame(metadata),by="sequenceID")
  ageClassI<-ageClassI[!duplicated(ageClassI$wormbaseID),]
  ageClassII<-ageClassII[!duplicated(ageClassII$wormbaseID),]
  write.table(ageClassI,paste0(outPath,"/publicData/agingClassI_Murphy2003.csv"),
              row.names=F,quote=T,sep=";")
  write.table(ageClassII,paste0(outPath,"/publicData/agingClassII_Murphy2003.csv"),
              row.names=F,quote=T,sep=";")
  file.remove(paste0(outPath,"/publicData/41586_2003_BFnature01789_MOESM2_ESM.xls"))
}



#######################-
## Aging txptome Tarkhov(2019)------
#######################-
# A universal transcriptomic signature of age reveals the temporal scaling of Caenorhabditis elegans aging trajectories
# Andrei E. Tarkhov, Ramani Alla, Srinivas Ayyadevara, Mikhail Pyatnitskiy, Leonid I. Menshikov, Robert J. Shmookler Reis & Peter O. Fedichev
# Scientific Reports volume 9, Article number: 7368 (2019)
# https://www.nature.com/articles/s41598-019-43075-z#Abs1
downSignatrURL<-"https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-019-43075-z/MediaObjects/41598_2019_43075_MOESM6_ESM.csv"
upSignatrURL<-"https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-019-43075-z/MediaObjects/41598_2019_43075_MOESM7_ESM.csv"
ageBiomrkrURL<-"https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-019-43075-z/MediaObjects/41598_2019_43075_MOESM8_ESM.csv"

if(!file.exists(paste0(outPath,"/publicData/agingBiomarkers_Tarkhov2019.csv"))){
  download.file(downSignatrURL,paste0(outPath,"/publicData/41598_2019_43075_MOESM6_ESM.csv"))
  download.file(upSignatrURL,paste0(outPath,"/publicData/41598_2019_43075_MOESM7_ESM.csv"))
  download.file(ageBiomrkrURL,paste0(outPath,"/publicData/41598_2019_43075_MOESM8_ESM.csv"))
  # there is some in the file. need to open manually and save as xlsx
  downSignatr<-read.csv(paste0(outPath,"/publicData/41598_2019_43075_MOESM6_ESM.csv"),
                        sep=",",header=T,skip=1)
  downSignatr$X<-NULL
  colnames(downSignatr)<-c("wormbaseID","description","p.value")
  upSignatr<-read.csv(paste0(outPath,"/publicData/41598_2019_43075_MOESM7_ESM.csv"),
                      sep=",",header=T,skip=1)
  upSignatr$X<-NULL
  colnames(upSignatr)<-c("wormbaseID","description","p.value")
  ageBiomrkr<-read.csv(paste0(outPath,"/publicData/41598_2019_43075_MOESM8_ESM.csv"),
                       sep=",",header=T,skip=1)
  ageBiomrkr$X<-NULL
  colnames(ageBiomrkr)<-c("wormbaseID","description","Coefficient")

  downSignatr<-inner_join(downSignatr,as.data.frame(metadata),by="wormbaseID") #67
  upSignatr<-inner_join(upSignatr,as.data.frame(metadata),by="wormbaseID") #260
  ageBiomrkr<-inner_join(ageBiomrkr,as.data.frame(metadata),by="wormbaseID") #71

  write.table(downSignatr,paste0(outPath,"/publicData/agingDownSignature_Tarkhov2019.csv"),
              row.names=F,quote=T,sep=";")
  write.table(upSignatr,paste0(outPath,"/publicData/agingUpSignature_Tarkhov2019.csv"),
              row.names=F,quote=T,sep=";")
  write.table(ageBiomrkr,paste0(outPath,"/publicData/agingBiomarkers_Tarkhov2019.csv"),
              row.names=F,quote=T,sep=";")

  file.remove(paste0(outPath,"/publicData/41598_2019_43075_MOESM6_ESM.csv"))
  file.remove(paste0(outPath,"/publicData/41598_2019_43075_MOESM7_ESM.csv"))
  file.remove(paste0(outPath,"/publicData/41598_2019_43075_MOESM8_ESM.csv"))
}


######################3-
## Aging Riedel (2013) (Murphy lab) -----
######################-
#https://www.nature.com/articles/ncb2720
#DAF-16 employs the chromatin remodeller SWI/SNF to promote stress resistance and longevity
# Christian G. Riedel, Robert H. Dowen, Guinevere F. Lourenco, Natalia V. Kirienko, Thomas Heimbucher, Jason A. West, Sarah K. Bowman, Robert E. Kingston, Andrew Dillin, John M. Asara & Gary Ruvkun
# Nature Cell Biology volume 15, pages 491–501 (2013)
# data not available in supplementary material... need to do alignment and deseq2

# daf-2
# daf-2/daf-16
# got table of metadata from SRA run selector:
df<-read.delim(paste0(outPath,"/publicData/SraRunTable_Riedel2013_SRP017908.txt"))
#df<-df[,c("Run","SRA_Study","Sample_Name")]
df$replicate<-stringr::str_sub(df$Sample_Name,stringr::str_locate(df$Sample_Name,"r.?$"))
df$Sample_Name<-gsub("_d1a_mRNA_r.?$","",df$Sample_Name)
df1<-df[,c("Run","SRA_Study","Sample_Name","replicate","strain")]
colnames(df1)<-c("SRRnumber","dataset","bioType","replicate","strain")
write.table(df1,file=paste0(outPath,"/publicData/SRR_Riedel2013_SRP017908.tsv"),row.names=F,quote=F)

######################3-
## Aging Zarse (2012) (Ristow lab) ------
######################-
# Cell Metabolism
# Volume 15, Issue 4, 4 April 2012, Pages 451-465
# Impaired Insulin/IGF1 Signaling Extends Life Span by Promoting Mitochondrial L-Proline Catabolism to Induce a Transient ROS Signal
# KimZarse12SebastianSchmeisser13MarcoGroth4SteffenPriebe5GregorBeuster1DoreenKuhlow16ReinhardGuthke5MatthiasPlatzer4C. RonaldKahn2MichaelRistow16
# https://www.sciencedirect.com/science/article/pii/S1550413112000940?via%3Dihub#app2

#daf-2
df<-read.csv(paste0(outPath,"/publicData/SraRunTable_Zarse2012_GSE36041.txt"))
df<-df[df$Organism=="Caenorhabditis elegans",]
df$replicate<-stringr::str_sub(df$Library.Name,stringr::str_locate(df$Library.Name,".$"))
df$Genotype<-gsub("wild type","wt",df$Genotype)
df$Genotype<-gsub("daf-2\\(e1370\\)","daf-2",df$Genotype)
df1<-df %>% dplyr::group_by(SRA.Study,Genotype,replicate) %>% dplyr::summarise(SRRnum=paste(Run,collapse=";"))
colnames(df1)<-c("dataset","bioType","replicate","SRRnumber")
df1<-df1[c(4,1,2,3)]
write.table(df1,file=paste0(outPath,"/publicData/SRR_Zarse2012_GSE36041.tsv"),row.names=F,quote=F)


######################-
## Aging Heestand (2013) (Antebi lab) ------
######################-
# Dietary Restriction Induced Longevity Is Mediated by Nuclear Receptor NHR-62 in Caenorhabditis elegans
# Bree N. Heestand, Yidong Shen, Wei Liu, Daniel B. Magner, Nadia Storm, Caroline Meharg, Bianca Habermann, Adam Antebi
# Published: July 25, 2013
# https://doi.org/10.1371/journal.pgen.1003651

eat2URL<-"https://doi.org/10.1371/journal.pgen.1003651.s012"
eat2File<-paste0(outPath,"/publicData/journal.pgen.1003651.s012.xlsx")
if(remakeFiles){
  file.remove(paste0(outPath,"/publicData/eat2down_Heestand2012.csv"))
}

if(!file.exists(paste0(outPath,"/publicData/eat2down_Heestand2012.csv"))){
  download.file(eat2URL,destfile= eat2File)
  eat2<-read_excel(paste0(outPath,"/publicData/journal.pgen.1003651.s012.xlsx"))
  eat2<-eat2[,c("Gene ID","eat-2 vs N2 fold-change (log2)","eat-2 vs N2 (padj)")]
  colnames(eat2)<-c("sequenceID","eat2_lfc","eat2_padj")
  eat2up<-getSignificantGenes(eat2,padj=0.05,lfc=2,namePadjCol="eat2_padj",
                              nameLfcCol="eat2_lfc",direction="gt",
                              chr="all") #361
  write.csv(eat2up,paste0(outPath,"/publicData/eat2up_Heestand2012.csv"),
            quote=F,row.names=F)

  eat2down<-getSignificantGenes(eat2,padj=0.05,lfc= -2,
                                namePadjCol="eat2_padj",
                                nameLfcCol="eat2_lfc",direction="lt",
                                chr="all") #537
  write.csv(eat2down,paste0(outPath,"/publicData/eat2down_Heestand2012.csv"),
            quote=F,row.names=F)
  file.remove(paste0(outPath,"/publicData/journal.pgen.1003651.s012.xlsx"))
}


######################-
## Aging Ayyadevara (2009) (Shmookler Reis lab) ------
######################-
# Caenorhabditis elegans PI3K mutants reveal novel genes underlying exceptional stress resistance and lifespan
# Srinivas Ayyadevara, Çagdaþ Tazearslan, Puneet Bharill, Ramani Alla, Eric Siegel, Robert J. Shmookler Reis,
# First published: 17 November 2009
#https://onlinelibrary.wiley.com/doi/10.1111/j.1474-9726.2009.00524.x



#McElwee (2003) paper
# McElwee, J., Bubb, K., and Thomas, J.H. (2003). Transcriptional outputs of the Caenorhabditis elegans forkhead protein DAF-16. Aging Cell 2, 111–121.
# https://www.deepdyve.com/lp/wiley/transcriptional-outputs-of-the-caenorhabditis-elegans-forkhead-protein-T0bliMYws0


######################-
## gene expression map (mountains, stuart kim)-----
######################-
# https://science.sciencemag.org/content/293/5537/2087/tab-figures-data




#####################-
## Broad expression vs regulated Gerstein et al. (2014) (modEncode)------
#####################-
#https://www.encodeproject.org/comparative/transcriptome/
wormGeneURL<-"http://cmptxn.gersteinlab.org/worm_gene.xlsx"

if(remakeFiles){
  file.remove(paste0(outPath,"/publicData/broadVregExpn_Gerstein2014.csv"))
}

if(!file.exists(paste0(outPath,"/publicData/broadVregExpn_Gerstein2014.csv"))) {
  download.file(wormGeneURL,paste0(outPath,"/publicData/worm_gene.xlsx"))
  wormGene<-readxl::read_excel(paste0(outPath,"/publicData/worm_gene.xlsx"))
  broad<-wormGene[,c("Gene","expr","BroadlyExpr_Score","L3_N2_L3-1")]
  colnames(broad)<-c("sequenceID","allExpr","broadExprScore","L3expr")
  hist(broad$broadExprScore,breaks=100)
  sum(!is.na(broad$broadExprScore)) # 5180 this is number of broadly expressed genes in supl therefore any gene htat has a score is broadly expressed. all otherse are regulated?
  #smoothScatter(log2(broad$L3expr),broad$broadExprScore,xlim=c(0,15),main="broad expression score vs L3 expression")
  #smoothScatter(log2(broad$allExpr),broad$broadExprScore,xlim=c(0,15),main="Broad expression score vs all stage average")
  ## not clear what threshold to use for expression
  broad$category<-NA
  broad$category[!is.na(broad$broadExprScore)]<-"broad"
  minBroadExpn<-min(broad$allExpr[broad$category=="broad"],na.rm=T)
  print(minBroadExpn)#1.559
  broad$category[broad$L3expr>minBroadExpn & is.na(broad$category)]<-"reg_L3"
  broad$category[broad$allExpr>minBroadExpn & is.na(broad$category)]<-"reg_nonL3"
  broad$category[is.na(broad$category)]<-"low"
  print(table(broad$category))
  #broadExpression   lowExpression    regulated_L3 regulated_nonL3
  #5180            6367            7990             840
  #    broad       low    reg_L3 reg_nonL3
  #   5180      7266      7120       811
  write.csv(broad,file=paste0(outPath,"/publicData/broadVregExpn_Gerstein2014.csv"),
            row.names=F,quote=F)
  file.remove(paste0(outPath,"/publicData/worm_gene.xlsx"))
}



#####################-
## Chromatin states Evans et al. (2016) (Ahinger lab)------
#####################-
# https://www.pnas.org/content/pnas/113/45/E7020.full.pdf

# Dataset S1. Coordinates of EE and L3 chromatin states
# File of chromosome, start position, end position, and state number. Coordinates are in WS220 and follow BED conventions (start positions are in zero-based coordinates and end positions in one-based coordinates).

if(remakeFiles | !file.exists(paste0(outPath,"/publicData/chromDomains_L3_Evans2016_ce11.bed"))){
  chromStatesURL<-"https://www.pnas.org/highwire/filestream/623778/field_highwire_adjunct_files/0/pnas.1608162113.sd01.xlsx"
  download.file(chromStatesURL,paste0(outPath,"/publicData/",basename(chromStatesURL)))
  ce10toCe11url<-"http://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz"
  ce10toCe11<-"ce10Toce11.over.chain"
  download.file(ce10toCe11url,paste0(outPath,"/publicData/",ce10toCe11,".gz"))
  system(paste0("gunzip ",outPath,"/publicData/",ce10toCe11,".gz"))
  file.remove(paste0(outPath,"/publicData/",ce10toCe11,".gz"))

  chrAstates<-readxl::read_excel(paste0(outPath,"/publicData/pnas.1608162113.sd01.xlsx"),sheet="L3 autosome states",col_names=c("chr","start","end","state"))
  chrXstates<-readxl::read_excel(paste0(outPath,"/publicData/pnas.1608162113.sd01.xlsx"),sheet="L3 chr X states",col_names=c("chr","start","end","state"))


  chrStates<-rbind(chrAstates,chrXstates)
  chrStatesGR<-GRanges(paste0("chr",chrStates$chr,":",
                              (chrStates$start+1),"-",
                              chrStates$end))
  chrStatesGR$score<-c(chrAstates$state,chrXstates$state)
  chainCe10toCe11<-import.chain(paste0(outPath,"/publicData/",ce10toCe11))
  chrStatesGR_ce11<-unlist(liftOver(chrStatesGR,chain=chainCe10toCe11))

  seqinfo(chrStatesGR_ce11)<-ce11seqinfo
  export(chrStatesGR_ce11,
         con=paste0(outPath,"/publicData/chromStates_L3_Evans2016_ce11.bed"),
         format="bed")
  file.remove(paste0(outPath,"/publicData/",basename(chromStatesURL)))

  # Dataset S2. Coordinates of EE and L3 domains
  # Excel file of chromosome, start position, end position of EE and L3 active domains, border regions, and regulated domains (each in separate tab). Additionally, border regions have strand information in column six to indicate if active domain is on the left (−) or on the right (+). Coordinates are in WS220 and follow BED conventions.
  chromDomainsURL<-"https://www.pnas.org/highwire/filestream/623778/field_highwire_adjunct_files/1/pnas.1608162113.sd02.xlsx"
  download.file(chromDomainsURL,paste0(outPath,"/publicData/",basename(chromDomainsURL)))

  l3active<-readxl::read_excel(paste0(outPath,"/publicData/",
                                      basename(chromDomainsURL)),
              sheet="L3 active domains",
              col_names=c("chr","start","end"))
  l3active$name<-"active"
  l3active$score<-"."
  l3active$strand<-"*"
  l3regulated<-readxl::read_excel(paste0(outPath,"/publicData/",
                                         basename(chromDomainsURL)),
                          sheet="L3 regulated domains",
                          col_names=c("chr","start","end"))
  l3regulated$name<-"regulated"
  l3regulated$score<-"."
  l3regulated$strand<-"*"
  l3borders<-readxl::read_excel(paste0(outPath,"/publicData/",basename(chromDomainsURL)),sheet="L3 borders",col_names=c("chr","start","end","name","score","strand"))
  l3borders$name<-"border"


  chrDomains<-rbind(l3active,l3regulated,l3borders)
  chrDomainsGR<-GRanges(paste0("chr",chrDomains$chr,":",(chrDomains$start+1),"-",
                               chrDomains$end,":",chrDomains$strand))
  mcols(chrDomainsGR)<-chrDomains[,c("name","score")]
  chainCe10toCe11<-import.chain(paste0(outPath,"/publicData/",ce10toCe11))
  chrDomainsGR_ce11<-unlist(liftOver(chrDomainsGR,chain=chainCe10toCe11))

  seqinfo(chrDomainsGR_ce11)<-ce11seqinfo
  chrDomainsGR_ce11$score<-NULL
  export(chrDomainsGR_ce11,
         con=paste0(outPath,"/publicData/chromDomains_L3_Evans2016_ce11.bed"),
         format="bed")
  file.remove(paste0(outPath,"/publicData/",basename(chromDomainsURL)))
}

###################-
## tissue specificity Kalezsky et al. 2018-----
###################-
# https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007559#pgen.1007559.s017
# Transcriptome analysis of adult Caenorhabditis elegans cells reveals tissue-specific gene and isoform expression
#Rachel Kaletsky , Victoria Yao , April Williams, Alexi M. Runnels, Alicja Tadych, Shiyi Zhou, Olga G. Troyanskaya , Coleen T. Murphy

### ubiquitously (?) expressed  - Table S1
ubiqURL<-"https://doi.org/10.1371/journal.pgen.1007559.s010"


### Tissue enriched  - Table S8
tissueEnrichURL<-"https://doi.org/10.1371/journal.pgen.1007559.s017"


### Tissue specific - Table S9
tissueSpecificURL<-"https://doi.org/10.1371/journal.pgen.1007559.s018"


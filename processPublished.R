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
library(ggVennDiagram)
library(ggplot2)
library(DESeq2)

source("./variableSettings.R")
source("./functions.R")

if(!dir.exists(paste0(outPath,"/publicData"))) {
  dir.create(paste0(outPath,"publicData"))
}

remakeFiles=TRUE
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


#######################
## manually curated from papers
#######################


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


#######################
## Jans 2009
#######################

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




#######################
## Kramer 2015
#######################
kramerURL<-"https://doi.org/10.1371/journal.pgen.1005698.s011"
kramerFileName="Kramer_2015_PlotGen_S3file.xlsx"
if(remakeFiles){
  file.remove(paste0(outPath,"/publicData/kramer2015_gr.rds"))
}

if(! file.exists(paste0(outPath,"/publicData/kramer2015_gr.rds"))){
  download.file(url=kramerURL,destfile=paste0(outPath,"/publicData/",kramerFileName))

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

  saveRDS(kramergr,file=paste0(outPath,"/publicData/kramer2015_gr.rds"))
  file.remove(paste0(outPath,"/publicData/",kramerFileName))
}

kramergr<-readRDS(file=paste0(outPath,"/publicData/kramer2015_gr.rds"))

lfcVal<-0.5
padjVal<-0.05
idx<-!is.na(kramergr$dpy27_RNAi_L3_padj) &
  kramergr$dpy27_RNAi_L3_padj < padjVal &
  kramergr$dpy27_RNAi_L3_log2_fold_change > lfcVal &
  seqnames(kramergr)=="chrX"
kramerdpy27dc<-kramergr[idx,]
colIdx<-grep("(set)|(dpy21)|(mixed_sex)",colnames(mcols(kramerdpy27dc)))
mcols(kramerdpy27dc)[,colIdx]<-NULL
saveRDS(kramerdpy27dc,file=paste0(outPath,"/publicData/kramer2015_chrXup_dpy27_lfc",
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
saveRDS(kramerdpy21dc,file=paste0(outPath,"/publicData/kramer2015_chrXup_dpy21_lfc",
                                  formatC(lfcVal,format="e",digits=0),"_p",
                                  formatC(padjVal,format="e",digits=0),
                                  "_gr.rds"))




###############-
#  Meeuse et al 2020 - Oscillating genes ----------------------------------
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
#  Latorre et al 2015 - Oscillating genes ----------------------------------
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

osc<-read.delim(paste0(outPath,"/publicData/oscillatingGenes.tsv"))
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

  hsUP<-getSignificantGenes(garrigues, padj=0.05, lfc=1,
                            namePadjCol="P-adj",
                            nameLfcCol="log2(FC)", direction="gt")
  hsUP
  saveRDS(hsUP,file=paste0(outPath,"/publicData/hsUp_garrigues2019.rds"))

  hsDOWN<-getSignificantGenes(garrigues, padj=0.05, lfc=-1,
                              namePadjCol="P-adj",
                              nameLfcCol="log2(FC)", direction="lt")
  hsDOWN
  saveRDS(hsDOWN,file=paste0(outPath,"/publicData/hsDown_garrigues2019.rds"))

  file.remove(paste0(outPath,"/publicData/",garriguesFileName))
}



###############################
## get germline-soma genes from Boeck-Waterston_GR2016
###############################

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

  germline<-resLFCgermline[resLFCgermline$padj<0.05 & resLFCgermline$log2FoldChange>0.5,]
  nongl<-resLFCgermline[resLFCgermline$padj<0.05 & resLFCgermline$log2FoldChange< -0.5,]
  germline$wormbaseID<-rownames(germline)
  nongl$wormbaseID<-rownames(nongl)
  dim(germline) #2749
  dim(nongl) #4993


  resLFCsoma<- lfcShrink(dds,contrast=c("germline","soma","mix"),type="ashr")
  plotMA(resLFCsoma,main="soma")

  resLFCsoma<-resLFCsoma[!is.na(resLFCsoma$padj),]
  soma<-resLFCsoma[resLFCsoma$padj<0.05 & resLFCsoma$log2FoldChange>0.5,]
  nonsoma<-resLFCsoma[resLFCsoma$padj<0.05 & resLFCsoma$log2FoldChange< -0.5,]
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


###############################
## get germline-soma genes from Reinke_DEV2004
###############################

glFile="germlineVsoma_Reinke_Dev2004"

if(remakeFiles){
  file.remove(paste0(outPath,"/publicData/germlineGenes_Reinke2004.csv"))
}

if (!file.exists(paste0(outPath,"/publicData/germlineGenes_Reinke2004.csv"))) {
  link2="http://dev.biologists.org/highwire/filestream/1201187/field_highwire_adjunct_files/0/Data_S1.zip"
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



######################
## compare germline datasets
######################

glData<-read.csv(paste0(outPath,"/publicData/germlineGenes_Reinke2004.csv"))
glvSoma<-read.csv(paste0(outPath,"/publicData/germlineSomaGenes_Boeck2016.csv"))
x<-list(germline=glData$wormbaseID,
        germlineL4=glvSoma$wormbaseID[glvSoma$germline=="germlineL4"],
        somaL4=glvSoma$wormbaseID[glvSoma$germline=="somaL4"])

p1<-ggVennDiagram(x) + ggtitle(label=paste0("Reinke(2004) vs Boeck(2016) L4 soma"))

x<-list(germline=glData$wormbaseID,
        gonadYA=glvSoma$wormbaseID[glvSoma$germline=="gonadYA"],
        somaYA=glvSoma$wormbaseID[glvSoma$germline=="somaYA"])

p2<-ggVennDiagram(x) + ggtitle(label=paste0("Reinke(2004) vs Boeck(2016) YA gonad"))

x<-list(germlineL4=glvSoma$wormbaseID[glvSoma$germline=="germlineL4"],
        gonadYA=glvSoma$wormbaseID[glvSoma$germline=="gonadYA"],
        somaL4=glvSoma$wormbaseID[glvSoma$germline=="somaL4"],
        somaYA=glvSoma$wormbaseID[glvSoma$germline=="somaYA"])

p3<-ggVennDiagram(x) + ggtitle(label=paste0("Boeck(2016) L4 vs YA soma/germline"))

p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
ggplot2::ggsave(filename=paste0(outPath, "/publicData/venn_ReinkeVsBoeck.pdf"),
                plot=p, device="pdf",width=29,height=11,units="cm")

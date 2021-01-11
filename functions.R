#' Make directories
#'
#' @param path String with path to where the directories should be made
#' @param dirNameList Vector of strings with names of directories to create (can include multilevel directories)
#' @return Creates the directories listed in dirNameList
#' @examples
#' makeDirs(path=".",dirNameList=c("/txt","/rds/sample1"))
#' @export
makeDirs<-function(path,dirNameList=c()) {
  sub("\\/$","",path) #remove directory slash if present in path string
  for (d in dirNameList) {
    if (!dir.exists(paste0(path,"/",d))){  # for alignments
      dir.create(paste0(path,"/",d), recursive=TRUE, showWarnings=FALSE)
    }
  }
}


#' Filter DESeq2 table results
#'
#' @param resultsTable Table of DESeq results
#' @param padj Adjusted p value threshold
#' @param lfc Log fold change threshold
#' @param direction Whether to find genes that are less than (lt), or greater than (gt) the log fold change threshold, or both extreme tails ("both")
#' @param chr Include all genes in genome ("all") only those on the X chromosome ("chrX"), or only autosomes ("autosomes")
#' @param outPath Path to working directory
#' @param filenamePrefix Text to add to filename
#' @param writeTable Should results table be automatically saved to a file (default: True)
#' @param IDcolumnName Name of column containing feature IDs (default: "ID")
#' @return filtered table of results which is also automatically written to disk
#' @export
filterResults<-function(resultsTable, padj=0.05, lfc=0, direction="both",
                        chr="all", outPath=".", filenamePrefix="",
                        writeTable=T,IDcolName="ID") {
  sigGenes<-getSignificantGenes(resultsTable,padj,lfc,direction=direction,chr=chr)
  idx<-resultsTable[,IDcolName] %in% sigGenes[,IDcolName]
  filtTable<-resultsTable[idx,c("baseMean","log2FoldChange","padj",
                                IDcolName,"chr","start","end","strand")]
  if(writeTable){
    if(!dir.exists(paste0(outPath,"/txt"))){
      dir.create(paste0(outPath,"/txt"))
    }
    write.csv(filtTable,file=paste0(outPath,"/txt/filtResults_p",
                                  padj,"_",direction,"-","lfc",
                                  lfc,"_",chr,".csv"), row.names=F,
            quote=F)
  }
  return(filtTable)
}




#' Get significant genes from  RNAseq results
#'
#' @param resultsTable Table of DESeq results
#' @param padj Adjusted p value threshold
#' @param lfc Log fold change threshold
#' @param namePadjCol Name of column with adjusted P values
#' @param nameLFCcol Name of column with log fold change values
#' @param direction Whether to find genes that are less than (lt), or greater than (gt) the log fold change threshold, or both extreme tails ("both")
#' @param chr Include all genes in genome ("all") only those on the X chromosome ("chrX"), or only autosomes ("autosomes")
#' @param nameChrCol Name of column with chromosome names.
#' @return Filtered table of significant genes at a certain log fold change and adjusted p value.
#' @export
getSignificantGenes<-function(resultsTable, padj=0.05, lfc=0, namePadjCol="padj",
                              nameLfcCol="log2FoldChange", direction="both",
                              chr="all", nameChrCol="chr", outPath="."){
  if(direction=="both") {
    idx<-!is.na(resultsTable[,namePadjCol]) & resultsTable[,namePadjCol]<padj & abs(resultsTable[,nameLfcCol])>lfc
  } else if(direction=="gt") {
    idx<-!is.na(resultsTable[,namePadjCol]) & resultsTable[,namePadjCol]<padj & resultsTable[,nameLfcCol]>lfc
  } else if(direction=="lt") {
    idx<-!is.na(resultsTable[,namePadjCol]) & resultsTable[,namePadjCol]<padj & resultsTable[,nameLfcCol]<lfc
  } else {
    print("direction must be 'both' to get both tails, \n'gt' to get lfc larger than a specific value, \nor 'lt' to get lfc less than a certain value")
  }
  if(chr=="all"){
    idx<-idx
  } else if(chr=="chrX"){
    idx<-idx & !is.na(resultsTable[,nameChrCol]) & resultsTable[,nameChrCol]=="chrX"
  } else if(chr=="autosomes"){
    idx<-idx & !is.na(resultsTable[,nameChrCol]) & resultsTable[,nameChrCol]!="chrX"
  } else {
    print("chr must be one of 'all', 'chrX' or 'autosomes'")
  }
  filtTable<-resultsTable[idx,]
  return(filtTable)
}


#' Convert short gene name to proper format
#'
#' Takes a gene name like scc1cs and adds hyphen. Function is based
#' on finding one or more digits in the name and putting a hyphen
#' in front of it.
#' @param geneName Gene name without hyphen, e.g. scc1cs
#' @return Gene name with hyphen, e.g. scc-1cs
#' @export
prettyGeneName<-function(geneName){
  gsub("([[:digit:]]+)","-\\1",geneName)
}


#' Assign GRanges to A/B compartment
#'
#' Function takes in a genomic ranges object which you wish to assign,
#' and a genomic ranges for an eigen vector giving the compartments,
#' and then adds pcaScore and compartment (A or B) columns to the
#' metadata of the GRanges object. If strings for the names of the gr
#' and pca are given a bedgraph will be produced with the AB
#' assignments of all the genes.
#' @param gr Genomic ranges object which you wish to assign to A/B compartments
#' @param pcagr Genomic ranges for eigen vector giving the A/B compartments
#' @param grName Name of the data in the gr object. Used to output a bedgraph of AB assignments
#' @param pcaName Name of the data in the pcagr object. Used to output a bedgraph of AB assignments
#' @return The gr genomic ranges that was input with additional metadata columns of pcaScore and compartment
assignGRtoAB<-function(gr, pcagr,grName=NULL,pcaName=NULL,
                       outPath="."){
  ol<-as.data.frame(findOverlaps(gr,pcagr,ignore.strand=T))
  ol$subjectScore<-pcagr$score[ol$subjectHits]
  pcaScore<-ol %>% group_by(queryHits)%>% dplyr::summarise(pcaScore=mean(subjectScore,na.rm=T))

  gr$pcaScore[pcaScore$queryHits]<-pcaScore$pcaScore
  gr$compartment<-as.factor(ifelse(gr$pcaScore>0,"B","A"))
  idx<-is.na(gr$pcaScore)
  print(paste0(sum(idx)," genes have no overlapping PCA bin"))
  gr<-gr[! idx ]
  forBG<-gr
  forBG$score<-ifelse(forBG$compartment=="A",1,-1)
  if(!(is.null(grName) | is.null(pcaName))){
    export(forBG,con=paste0(outPath,"/tracks/",grName,"_",
                            "__Compartments_",pcaName, ".bedGraph"),
           format="bedGraph")
  }
  return(gr)
}

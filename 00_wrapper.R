### make sure to put the correct settings in variableSettings.R

source("variableSettings.R")

# make metadata object
if(!file.exists(paste0(outPath,"/wbGeneGR_WS275.rds"))){
  source("./createMetadataObj.R")
}

# process some public datasets required later (only do once)
print("processing published datasets")
source("processPublished.R")

# run core DESeq2 pipeline
print("main DESeq2 analysis")
if(combineChrAX){
  source("combine_chrAchrX.R")
} else {
  source("DESeqAnalysis_SALMON.R")
}

#
#make bigwigs from STAR data
#source("makeSTARbw.R")

# find best LFC threshold
print("ROCit plots for optimal LFC threshold")
source("ROCit.R")

# get read counts for each sample
print("collecting read counts from qc files")
source("collectAllCountData.R")

# compare different biological interventions to eachother
print("comparing the datasets to one another")
source("compareDatasets.R")

# compare pval and lfc for specific lists of genes
print("Check results for particular lists of genes")
source("compareGeneLists.R")

# compare HiC features - domains, loops
print("compare to HiC features")
source("compareHICfeatures.R")

# compare tissue and GO enrichment
print("check tissue and GO enrichment")
source("compareTissueAndGO.R")

if(filterData){
  fileNamePrefix<-filterPrefix
}
rmarkdown::render(input="AllPlots_TEA.Rmd",output_format="pdf_document",
                  output_file=paste0(outPath,"/tissue/tea/",fileNamePrefix,"allPlots_TEA.pdf"))

rmarkdown::render(input="AllPlots_WORMCAT.Rmd", output_format="pdf_document",
                  output_file=paste0(outPath,"/wormcat/",fileNamePrefix,"allPlots_WORMCAT.pdf"))

# compare KEGG pathway enrichment
print("check KEGG pathways")
source("compareKEGG.R")

# compare to dosage compensation data sets
print("compare lengths of regulated genes")
source("compareGeneLengths.R")

# compare to dosage compensation data sets
print("compare to public DC data sets")
source("compareToDCdatasets.R")

# check percentage germline vs soma genes among positives
print("look for soma/germline enrichment")
source("compareGermline.R")

# check enrichment for aging genes
print("look for longevity/aging enrichment")
source("compareAging.R")


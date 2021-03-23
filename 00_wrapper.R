### make sure to put the correct settings in variableSettings.R

# process some public datasets required later (only do once)
print("processing published datasets")
source("processPublished.R")
#
#make bigwigs from STAR data
#source("makeSTARbw.R")

# get read counts for each sample
#print("collecting read counts from qc files")
#source("collectAllCountData.R")

# run core DESeq2 pipeline
print("main DESeq2 analysis")
if(combineChrAX){
  source("combine_chrAchrX.R")
} else {
  source("DESeqAnalysis_SALMON.R")
}

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

# compare to dosage compensation data sets
print("compare to public DC data sets")
source("compareToDCdatasets.R")

# check percentage germline vs soma genes among positives
print("look for soma/germline enrichment")
source("compareGermline.R")


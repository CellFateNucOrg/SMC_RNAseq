
#source("processPublished.R")

# get read counts for each sample
source("collectAllCountData.R")

# run core DESeq2 pipeline
source("DESeqAnalysis_SALMON.R")

# compare different biological interventions to eachother
source("compareDatasets.R")

# compare pval and lfc for specific lists of genes
source("compareGeneLists.R")

# compare HiC features - domains, loops
source("compareHICfeatures.R")

# compare tissue and GO enrichment
source("compareTissueAndGO.R")

# compare to dosage compensation data sets
source("compareToDCdatasets.R")


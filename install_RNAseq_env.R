install.packages(c("ggplot2", "RColorBrewer","PoiClaClu","pheatmap","tidyr",
                 "gplots","ggpubr"))
install.packages(c("reticulate","xlsx","magick","cowplot"))

BiocManager::install(c("DESeq2","AnnotationDbi","Organism.dplyr","GenomicRanges",
                     "BSgenome.Celegans.UCSC.ce11", "tximport",
                     "GenomicFeatures","affy","EnhancedVolcano","apeglm",
                     "ashr","genomation"))


# for tissue enrichment
devtools::install_github("trinker/plotflow")
devtools::install_github("jsemple19/Wormcat")

reticulate::conda_install(envname="tea",packages="tissue_enrichment_analysis",pip=T, pip_options="git+https://github.com/dangeles/TissueEnrichmentAnalysis.git")
## note: there is a bug in this version of tea which need to be corrected in the
## code. Go to this file: vi ~/miniconda3/envs/tea/bin/tea and on line 66 change args.tisse_dictionary to args.dictionary





# Dpy26csRNAseq

Code for processing RNAseq data from Julie with strains 493 (dpy-26::TEVcs;hs::TEV) vs 500 (hs::TEV).

## Pipeline

On the server the data is mapped with the script:

_**mapRNAreads.sh3**_ using STAR and Salmon. 

The output is then analysed in R with DESeq2 using the following script:

**_DESeqAnalysis_JulieData_SALMON.R_**

## Comparison with other aligners

For comparison I also run similar scripts with the STAR output, or with QuasR:

**_DESeqAnalysis_JulieData_STAR.R_**

**_DESeqAnalysis_JulieData_QUASR.R_**

and then compared the output of the different methods using:

_**compareAligners.R**_

## Comparison with other datasets

Finally I compared the output of Salmon with other publically available datasets of genes with:

_**compareDatasets.R**_



## Publically available datasets

These datasets were retrieved manually:

**DCgenes_published.xlsx** 

**Jans2009_DC_suplTable4.txt**

**Jans2009_notDC_suplTable5.txt** 

**Kramer et al. (2015)** was retrieved automatcally.

All were processed with:

_**processPublished.R**_

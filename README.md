# SMC_RNAseq

Code for processing RNAseq data from Moushumi with strains wPM382 (dpy-26::TEVcs;hs::TEV), wPM784 (scc-1::TEVcs;hs::TEV) and wPM775 (kle-2::TEVcs;hs::TEV) vs 366 (hs::TEV).

The pipeline performs genomic alignment with STAR so that we can get bigwig files, and transcriptome alignment with Salmon in order to get better gene-level counts. Signficance is evaluated with DESeq2.

## Pipeline

When running the pipeline for the first time (use ubelix), the _**indexGenomeTranscripts.sh**_ needs to be run to index both the genome for STAR and the transcriptome for salmon.

On the server the data is mapped with the script:

_**mapRNAreads.sh**_ using STAR and Salmon. 

This requires a file **fastqList.txt** (adapt the **fastqList_example.txt** file) with three columns: fileName (full path), sampleName and repeatNum. processed files will be named using the sampleName and repeatNum.

The output is then analysed in R with DESeq2 using the following script:

**_DESeqAnalysis_SALMON.R_**



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

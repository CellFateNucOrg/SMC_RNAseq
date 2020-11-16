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

This produces the following plots:

- basic QC plots: bar plots and density plots of the raw counts,sample-sample clustering heatmap, heatmap of top 500 expressed genes with and without sample clustering, PCA plots coloured by main variables (_salmon_sampleQC.pdf_)

- boxplots of lfc of genes by chr or chr type (_salmon_xxxx_boxplots_expnByChr.pdf_, _salmon_xxxx_boxplots_expnByChrType.pdf_)

- MAplots (_salmon_xxxx_MAplots_results.pdf_)

- Volcano plots (_salmon_xxxxx_volcanoPlot_xxxxx.pdf_)

- heirarchical clustering of samples by most changed genes (_salmon_xxxxxx_hclust_mostChanged.pdf_)

- barplots of expression of a few individual genes, most changed in each sample (_salmon_xxxx_topGenes_normCounts.pdf_)

## Comparison of the different data sets

The significant genes from each data set are compared by:
_**compareDatasets.R**_

This script produces the following plots:

- barplots to campare counts of genes per chromosome (_./plots/bar_countsPerChr_xxxxxx.pdf_)

- venn diagrams to compare overlap between datasets (_./plots/venn_xxxxxx.pdf_)

- correlation of the log2 fold change of genes shared between data sets (_./plots/cor_xxxxx.png_)

## Comparison of HiC features

The association of gene expression with different HiC features are compared by:
_**compareHICfreatures.R**_


## Comparison of Tissue and GO term enrichment

_**compareTissueAndGO.R**_


## Publically available datasets

These datasets were retrieved manually:

**DCgenes_published.xlsx** 

**Jans2009_DC_suplTable4.txt**

**Jans2009_notDC_suplTable5.txt** 

**Kramer et al. (2015)** was retrieved automatcally.

All were processed with:

_**processPublished.R**_

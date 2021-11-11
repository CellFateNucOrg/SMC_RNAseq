#!/usr/bin/env bash

conda create --name RNAseq python=3.8
# OR in a specific location
#conda create --prefix /Volumes/Seagate/opt/miniconda3/env/salmon python=3.8

#conda activate ceftall
source $CONDA_ACTIVATE RNAseq

# OR in a specific location
#conda activate /Volumes/Seagate/opt/miniconda3/env/salmon
conda install cutadapt fastqc salmon star samtools htseq wiggletools -c bioconda -c conda-forge
conda install -c iuc ucsc_tools

conda create --name cufflinks python=3.5
source $CONDA_ACTIVATE cufflinks
conda install cufflinks -c bioconda

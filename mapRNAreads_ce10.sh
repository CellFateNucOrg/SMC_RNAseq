#! /usr/bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="RNAseq"
#SBATCH --time=0-8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=all
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-4

module add vital-it
module add UHTS/Quality_control/fastqc/0.11.7;
module add UHTS/Quality_control/cutadapt/2.5;
module add UHTS/Aligner/STAR/2.7.3a;
module add UHTS/Analysis/samtools/1.10;
export SALMON_SING="singularity exec /software/singularity/containers/salmon-1.2.1-1.ubuntu18.sif"
module add R/3.6.1;
module add UHTS/Aligner/bwa/0.7.17;
#module add UHTS/Analysis/HTSeq/0.9.1;
source $CONDA_ACTIVATE htseq

echo "current date"
date

echo "current git branch"
git branch

echo "current git version"
git log -1 --format="%H"

fastqFileList=./fastqList.txt
fastqFile=(`cut -f1 $fastqFileList`)
sampleNames=(`cut -f2 $fastqFileList`)
repeatNums=(`cut -f3 $fastqFileList`)
#laneNums=(`cut -f4 $fastqFileList`)

i=${SLURM_ARRAY_TASK_ID}

#####################
# mapping RNAseq data
#####################

###############################
########### VARIABLES #########
###############################

mRNAonly=false #false or true
fastqFile=${fastqFile[$i]}
sampleName=${sampleNames[$i]}
repeatNum=${repeatNums[$i]}
#laneNum=${laneNums[$i]}

nThreads=${SLURM_CPUS_PER_TASK}

genomeVer=WS220
GENOME_DIR=${HOME}/genomeVer/${genomeVer}
genomeFile=${GENOME_DIR}/sequence/c_elegans.${genomeVer}.genomic.fa
chromSizesFile=$GENOME_DIR/annotation/ws220.chrom.sizes


WORK_DIR=$PWD
QC_DIR=${WORK_DIR}/qc
WIGtoBW_DIR=${HOME}/mySoftware

#FASTQ_DIR=`dirname ${fastqFile}`
baseName=${sampleName}_${repeatNum}
#_${laneNum}

########################################################
### get initial read stats                            ##
########################################################
#
##run fastqc on sequences
#mkdir -p ${WORK_DIR}/qc/rawData
#fastqc ${fastqFile} -t $nThreads -o ${WORK_DIR}/qc/rawData 
#	
########################################################
### trim adaptors with cutadapt                       ##
########################################################
#
## use cutadapt to trim
#mkdir -p cutadapt
#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o cutadapt/${baseName}.fastq.gz -j $nThreads ${fastqFile}
#
##redo fastQC on trimmed reads
#mkdir -p ${WORK_DIR}/qc/cutadapt
#fastqc cutadapt/${baseName}.fastq.gz -t $nThreads -o ${WORK_DIR}/qc/cutadapt

########
# Align to genome with bwa aln
########
mkdir -p ${WORK_DIR}/bamBWA
echo "aligning $baseName to genome with BWA aln..."
bwa aln -t $nThreads $genomeFile ${WORK_DIR}/cutadapt/${baseName}.fastq.gz > ${WORK_DIR}/bamBWA/${baseName}.sai
bwa samse $genomeFile ${WORK_DIR}/bamBWA/${baseName}.sai ${WORK_DIR}/cutadapt/${baseName}.fastq.gz > ${WORK_DIR}/bamBWA/${baseName}.sam

echo "Removing unmapped reads and converting $baseName SAM to BAM..."
samtools view -b -F 4 -@ $nThreads ${WORK_DIR}/bamBWA/${baseName}.sam > ${WORK_DIR}/bamBWA/${baseName}.bam

rm ${WORK_DIR}/bamBWA/${baseName}.sam
rm ${WORK_DIR}/bamBWA/${baseName}.sai


#gtfFile=${GENOME_DIR}/annotation/c_elegans.PRJNA13758.${genomeVer}.annotations_rpt.gtf
annotFile_rpt=${GENOME_DIR}/annotation/c_elegans.${genomeVer}.annotations_rpt.gtf

# count reads per feature with htseq
mkdir -p ${WORK_DIR}/htseq
htseq-count -f bam -r name -a 0 -m union --nonunique random  ${WORK_DIR}/bamBWA/${baseName}.bam $annotFile_rpt > ${WORK_DIR}/htseq/${baseName}_union_random.txt
htseq-count -f bam -r name -a 0 -m union --nonunique none  ${WORK_DIR}/bamBWA/${baseName}.bam $annotFile_rpt > ${WORK_DIR}/htseq/${baseName}_union_none.txt
#--additional-attr=gene_name

samtools sort -T ${WORK_DIR}/bamBWA/${baseName}  -@ $nThreads -o ${WORK_DIR}/bamBWA/${baseName}_sort.bam  ${WORK_DIR}/bamBWA/${baseName}.bam
samtools index -@ $nThreads ${WORK_DIR}/bamBWA/${baseName}_sort.bam
rm ${WORK_DIR}/bamBWA/${baseName}.bam


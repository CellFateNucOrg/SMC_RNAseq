#! /usr/bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --partition=epyc2
#SBATCH --qos=job_epyc2
#SBATCH --time=0-04:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-8
#SBATCH --job-name="getSRR"
#SBATCH --cpus-per-task=2
##SBATCH --tmp=2G

######### Make sure to set number of array jobs to number of SRR numbers #######

#source ./envSettings.sh
# get working directory
#WORK_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
WORK_DIR=`pwd`
echo $WORK_DIR
FASTQ_DIR=$WORK_DIR/fastq
QC_DIR=$WORK_DIR/qc

# name of SRR file. Must have no header and four  tab separated columns.
# First column SRR number. Second column is the dataset name. Third column is biological type. Fourth column is replicate number
#SRRfile=${WORK_DIR}/SRR_McMurchy2017_all.tsv
SRRfile=${WORK_DIR}/publicData/SRR_Riedel2013_SRP017908.tsv

module load vital-it
module add UHTS/Quality_control/fastqc/0.11.7;
module add UHTS/Analysis/sratoolkit/2.10.7;

# SRRfile is a tab separated file with at least 4 columns: 
# First column has one or more SRR numbers of a single sample. when there are multiple SRR numbers they are separated by a semicolon
# Second column has the dataset name from GEO.
# Third column has a sample name that is common to all rows with the same biological type. i
# Fourth column contains a number for each replicate
# Fifth column contains some other variable that needs to be taken into account like batch, or condition

SRRfields=(`cut -f1 ${SRRfile}`)
datasets=(`cut -f2 ${SRRfile}`)
biotypes=(`cut -f3 ${SRRfile}`)
replicates=(`cut -f4 ${SRRfile}`)
#conditions=(`cut -f5 ${SRRfile}`)

mkdir -p $FASTQ_DIR
mkdir -p $QC_DIR

# there should be as many arrays tasks as there are SRRs
i=$SLURM_ARRAY_TASK_ID

#get sample + repeat number base name
baseName=${datasets[$i]}_${biotypes[$i]}_${replicates[$i]}   #_${conditions[$i]}
mkdir -p ${WORK_DIR}/${baseName}


# split SRRfields with multiple SRR numbers
IFS=';' read -r -a SRRs <<< "${SRRfields[$i]}"

# download data
for SRR in "${SRRs[@]}"
do
  echo ${baseNames[$i]}
  echo $SRR
  #prefetch -O ${WORK_DIR}/${baseName} -X 4G -o ${WORK_DIR}/${baseName}/${SRR}.srr  $SRR
  fasterq-dump -O ${WORK_DIR}/${baseName} -b 4G  -e $SLURM_CPUS_PER_TASK ${SRR}
done

# check if reads are paired end or not
#if [ -f "${WORK_DIR}/${SRR}.srr_1.fastq" ]; 
if [ -f "${WORK_DIR}/${SRR}_1.fastq" ];
then 
  isPE=true; 
else
  isPE=false;
fi


# merge multiple SRRs from same sample 
if ( "$isPE" )
then
  #cat ${WORK_DIR}/${baseName}/*.srr_1.fastq > ${FASTQ_DIR}/${baseName}_R1.fastq
  #cat ${WORK_DIR}/${baseName}/*.srr_2.fastq > ${FASTQ_DIR}/${baseName}_R2.fastq
  cat ${WORK_DIR}/${baseName}/*_1.fastq > ${FASTQ_DIR}/${baseName}_R1.fastq
  cat ${WORK_DIR}/${baseName}/*_2.fastq > ${FASTQ_DIR}/${baseName}_R2.fastq
  gzip -f ${FASTQ_DIR}/${baseName}_R1.fastq
  gzip -f ${FASTQ_DIR}/${baseName}_R2.fastq
else 
  #cat ${WORK_DIR}/${baseName}/*.srr.fastq > ${FASTQ_DIR}/${baseName}.fastq 
  cat ${WORK_DIR}/${baseName}/*.fastq > ${FASTQ_DIR}/${baseName}.fastq  
  gzip -f ${FASTQ_DIR}/${baseName}.fastq 
fi

#rm -rf ${WORK_DIR}/${baseName}

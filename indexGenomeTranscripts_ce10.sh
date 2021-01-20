#! /usr/bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="index_gentrome"
#SBATCH --time=0-04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=all
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1

module add vital-it
#module add UHTS/Aligner/STAR/2.6.0c;
module add UHTS/Aligner/STAR/2.7.3a;
#module add UHTS/Analysis/salmon/0.11.2;
export SALMON_SING="singularity exec /software/singularity/containers/salmon-1.2.1-1.ubuntu18.sif"
module add UHTS/Assembler/cufflinks/2.2.1
module add R/3.6.1
module add UHTS/Aligner/bwa/0.7.17;

nThreads=$SLURM_CPUS_PER_TASK

############
# prepare directories
############

genomeVer=WS220
GENOME_DIR=${HOME}/genomeVer/${genomeVer}
mkdir -p ${GENOME_DIR}/sequence
mkdir -p ${GENOME_DIR}/annotation

genomeFile=${GENOME_DIR}/sequence/c_elegans.${genomeVer}.genomic.fa
curl ftp://ftp.wormbase.org/pub/wormbase/releases/${genomeVer}/species/c_elegans/c_elegans.${genomeVer}.genomic.fa.gz -o ${genomeFile}.gz
gunzip ${genomeFile}.gz 

annotFile=${GENOME_DIR}/annotation/c_elegans.${genomeVer}.annotations.gff3.gz
curl ftp://ftp.wormbase.org/pub/wormbase/releases/${genomeVer}/species/c_elegans/c_elegans.${genomeVer}.annotations.gff3.gz -o $annotFile 
gunzip $annotFile

annotFile=${annotFile%.gz}
# use gffread from cufflinks to convert gff to gtf
b=(`basename -s .gff3 ${annotFile}`)
gffread $annotFile -T -o ${annotFile%/*}/${b}.gtf
#rm $annotFile

# need to remove wierd exons:
annotFile=${GENOME_DIR}/annotation/c_elegans.${genomeVer}.annotations.gtf
grep "exon" ${annotFile} > ${GENOME_DIR}/annotation/c_elegans.${genomeVer}.annotations1.gtf
mv  ${GENOME_DIR}/annotation/c_elegans.${genomeVer}.annotations1.gtf ${annotFile}
 
#mv c_elegans.PRJNA13758.${genomeVer}.annotations.gtf $annotFile


#####  index for repeats
# get repeat annotation from dfam
Rscript getRepeatData_ce10.R ${GENOME_DIR}/annotation

annotFile=${GENOME_DIR}/annotation/c_elegans.${genomeVer}.annotations.gtf
annotFile_rpt=${GENOME_DIR}/annotation/c_elegans.${genomeVer}.annotations_rpt.gtf
genomeFile=${GENOME_DIR}/sequence/c_elegans.${genomeVer}.genomic.fa
rptFile=${GENOME_DIR}/annotation/repeats_ce11_Dfam_2.0_nr.gtf

# combine coding transcriptome and repeat files
cat $annotFile $rptFile > $annotFile_rpt 


##### BWA index for repeats

# index genome
if [[ ! -f "${genomeFile}.bwt" ]]
then
  echo "indexing genome for bwa..."
  bwa index $genomeFile 
fi


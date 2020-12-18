#! /usr/bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="index_gentrome"
#SBATCH --time=0-08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=8G
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

mRNAonly=false
genomeVer=WS275
GENOME_DIR=${HOME}/genomeVer/${genomeVer}
mkdir -p ${GENOME_DIR}/sequence
mkdir -p ${GENOME_DIR}/annotation


#############
## index for STAR
#############
#
#genomeFile=${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.genomic.fa
#curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/genomic/c_elegans.PRJNA13758.${genomeVer}.genomic.fa.gz -o ${genomeFile}.gz
#gunzip ${genomeFile}.gz 
#
#annotFile=${GENOME_DIR}/annotation/c_elegans.PRJNA13758.${genomeVer}.annotations.gff3.gz
#curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/gff/c_elegans.PRJNA13758.${genomeVer}.annotations.gff3.gz -o $annotFile 
#gunzip $annotFile
#
#annotFile=${annotFile%.gz}
## use gffread from cufflinks to convert gff to gtf
#b=(`basename -s .gff3 ${annotFile}`)
#gffread $annotFile -T -o ${annotFile%/*}/${b}.gtf
##rm $annotFile
#
## need to remove wierd exons:
#annotFile=${GENOME_DIR}/annotation/c_elegans.PRJNA13758.${genomeVer}.annotations.gtf
#grep WormBase ${annotFile} > ${GENOME_DIR}/annotation/c_elegans.PRJNA13758.${genomeVer}.annotations1.gtf
#mv  ${GENOME_DIR}/annotation/c_elegans.PRJNA13758.${genomeVer}.annotations1.gtf ${annotFile}
# 
##mv c_elegans.PRJNA13758.${genomeVer}.annotations.gtf $annotFile
#
##index genome
##echo "indexing genome..."
#STAR --runMode genomeGenerate --genomeDir ${GENOME_DIR}/sequence --genomeFastaFiles ${genomeFile} --sjdbGTFfile ${annotFile} --runThreadN $nThreads --genomeSAindexNbases 12 
#
###### STAR index for repeats
## get repeat annotation from dfam
Rscript getRepeatData.R ${GENOME_DIR}/annotation

annotFile=${GENOME_DIR}/annotation/c_elegans.PRJNA13758.${genomeVer}.annotations.gtf
annotFile_rpt=${GENOME_DIR}/annotation/c_elegans.PRJNA13758.${genomeVer}.annotations_rpt.gtf
genomeFile=${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.genomic.fa
rptFile=${GENOME_DIR}/annotation/repeats_ce11_dfam_nr.gtf

cat $annotFile $rptFile > $annotFile_rpt 

#index genome
echo "indexing genome..."
mkdir -p ${GENOME_DIR}/sequence/repeats
STAR --runMode genomeGenerate --genomeDir ${GENOME_DIR}/sequence/repeats --genomeFastaFiles ${genomeFile} --sjdbGTFfile ${annotFile_rpt} --runThreadN $nThreads --genomeSAindexNbases 12



#############3
## index for Salmon
##############
#
#curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/genomic/c_elegans.PRJNA13758.${genomeVer}.genomic_masked.fa.gz -o ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.genomic_masked.fa.gz
#
#grep "^>" <(gunzip -c  ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.genomic_masked.fa.gz) | cut -d " " -f 1 > ${GENOME_DIR}/sequence/decoys.txt
#sed -i.bak -e 's/>//g' ${GENOME_DIR}/sequence/decoys.txt
#
#curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.${genomeVer}.mRNA_transcripts.fa.gz -o ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.mRNA_transcripts.fa.gz
#
#cat  ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.mRNA_transcripts.fa.gz ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.genomic_masked.fa.gz > ${GENOME_DIR}/sequence/mRNA_gentrome.fa.gz
#
#
#curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.${genomeVer}.ncRNA_transcripts.fa.gz -o ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.ncRNA_transcripts.fa.gz
#
#cat  ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.ncRNA_transcripts.fa.gz ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.genomic_masked.fa.gz > ${GENOME_DIR}/sequence/ncRNA_gentrome.fa.gz
#
#
#curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.${genomeVer}.pseudogenic_transcripts.fa.gz -o ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.pseudogenic_transcripts.fa.gz
#
#cat  ${GENOME_DIR}/sequence/c_elegans.PRJNA13758${genomeVer}.pseudogenic_transcripts.fa.gz ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.genomic_masked.fa.gz > ${GENOME_DIR}/sequence/pseudogenic_gentrome.fa.gz
#
#
#curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.${genomeVer}.transposon_transcripts.fa.gz -o ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.transposon_transcripts.fa.gz
#
#cat  ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.transposon_transcripts.fa.gz ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.genomic_masked.fa.gz > ${GENOME_DIR}/sequence/transposon_gentrome.fa.gz
#
##cd $GENOME_DIR
#
#${SALMON_SING} salmon index -t ${GENOME_DIR}/sequence/mRNA_gentrome.fa.gz -d ${GENOME_DIR}/sequence/decoys.txt -p $nThreads -i ${GENOME_DIR}/sequence/${genomeVer}_mRNA_index
#
#${SALMON_SING} salmon index -t ${GENOME_DIR}/sequence/ncRNA_gentrome.fa.gz -d ${GENOME_DIR}/sequence/decoys.txt -p $nThreads -i ${GENOME_DIR}/sequence/${genomeVer}_ncRNA_index
#
#${SALMON_SING} salmon index -t ${GENOME_DIR}/sequence/pseudogenic_gentrome.fa.gz -d ${GENOME_DIR}/sequence/decoys.txt -p $nThreads  -i ${GENOME_DIR}/sequence/${genomeVer}_pseudogenic_index
#
#${SALMON_SING} salmon index -t ${GENOME_DIR}/sequence/transposon_gentrome.fa.gz -d ${GENOME_DIR}/sequence/decoys.txt -p $nThreads -i ${GENOME_DIR}/sequence/${genomeVer}_transposon_index
#
## remove gentrome files
#rm ${GENOME_DIR}/sequence/mRNA_gentrome.fa.gz
#rm ${GENOME_DIR}/sequence/ncRNA_gentrome.fa.gz
#rm ${GENOME_DIR}/sequence/pseudogenic_gentrome.fa.gz
#rm ${GENOME_DIR}/sequence/transposon_gentrome.fa.gz
#
## use these addresses to reference indeces in main script
#mRNAindex=${GENOME_DIR}/sequence/${genomeVer}_mRNA_index
#ncRNAindex=${GENOME_DIR}/sequence/${genomeVer}_ncRNA_index
#pseudoIndex=${GENOME_DIR}/sequence/${genomeVer}_pseudogenic_index
#tnIndex=${GENOME_DIR}/sequence/${genomeVer}_transposon_index
#
#
#
#
################
### Index repeats
################
#
#
######## Salmon index for repeats
#curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.${genomeVer}.mRNA_transcripts.fa.gz -o ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.mRNA_transcripts.fa.gz
#
#grep "^>" <(gunzip -c  ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.mRNA_transcripts.fa.gz) | cut -d " " -f 1 > ${GENOME_DIR}/sequence/decoys_forRepeats.txt
#sed -i.bak -e 's/>//g' ${GENOME_DIR}/sequence/decoys_forRepeats.txt
#
#cat  ${GENOME_DIR}/annotation/repeats_ce11_dfam_nr.fa.gz ${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.mRNA_transcripts.fa.gz > ${GENOME_DIR}/sequence/repeats_gentrome.fa.gz
#
#kmerSize=15
#${SALMON_SING} salmon index -t ${GENOME_DIR}/sequence/repeats_gentrome.fa.gz -d ${GENOME_DIR}/sequence/decoys_forRepeats.txt -p $nThreads -i ${GENOME_DIR}/sequence/${genomeVer}_repeats_index_${kmerSize} --keepDuplicates -k $kmerSize
#rm ${GENOME_DIR}/sequence/repeats_gentrome.fa.gz
#
#repeatsIndex=${GENOME_DIR}/sequence/${genomeVer}_repeats_index_${kmerSize}



##### BWA index for repeats
genomeFile=${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.genomic.fa

# NOTE: cannot make gtf file for non-coding transcripts becuase of stringent gtf file format requirements.
# So will do star alignment to normal transcript index.
# will also try bwa alignment
# index genome
if [[ ! -f "${genomeFile}.bwt" ]]
then
  echo "indexing genome for bwa..."
  bwa index $genomeFile 
fi


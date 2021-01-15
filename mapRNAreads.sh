#! /usr/bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="RNAseq"
#SBATCH --time=0-8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=all
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1

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

genomeVer=WS275
GENOME_DIR=${HOME}/genomeVer/${genomeVer}
genomeFile=${GENOME_DIR}/sequence/c_elegans.PRJNA13758.${genomeVer}.genomic.fa
chromSizesFile=$GENOME_DIR/annotation/ws235.chrom.sizes
# annotFile=/home/ubelix/izb/semple/genomeVer/${genomeVer}/annotation/c_elegans.PRJNA13758.${genomeVer}.annotations.gff3.gz
# gunzip $annotFile
# annotFile=${annotFile%.gz}
# use gffread from cufflinks to convert gff to gtf
# module add UHTS/Assembler/cufflinks/2.2.1
# b=(`basename -s .gff3 ${annotFile}`)
# gffread $annotFile -T -o ${annotFile%/*}/${b}.gtf
# need to remove wierd exons:
# grep WormBase c_elegans.PRJNA13758.${genomeVer}.annotations.gtf > c_elegans.PRJNA13758.WS260.annotations1.gtf
# mv c_elegans.PRJNA13758.${genomeVer}.annotations1.gtf c_elegans.PRJNA13758.${genomeVer}.annotations.gtf
#mRNAseqFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.mRNA_transcripts.fa.gz
mRNAindex=${GENOME_DIR}/sequence/${genomeVer}_mRNA_index
ncRNAindex=${GENOME_DIR}/sequence/${genomeVer}_ncRNA_index
pseudoIndex=${GENOME_DIR}/sequence/${genomeVer}_pseudogenic_index
tnIndex=${GENOME_DIR}/sequence/${genomeVer}_transposon_index


WORK_DIR=$PWD
QC_DIR=${WORK_DIR}/qc
WIGtoBW_DIR=${HOME}/mySoftware

#FASTQ_DIR=`dirname ${fastqFile}`
baseName=${sampleName}_${repeatNum}
#_${laneNum}

#######################################################
## get initial read stats                            ##
#######################################################

#run fastqc on sequences
mkdir -p ${WORK_DIR}/qc/rawData
fastqc ${fastqFile} -t $nThreads -o ${WORK_DIR}/qc/rawData 
	
#######################################################
## trim adaptors with cutadapt                       ##
#######################################################

# use cutadapt to trim
mkdir -p cutadapt
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o cutadapt/${baseName}.fastq.gz -j $nThreads ${fastqFile}

#redo fastQC on trimmed reads
mkdir -p ${WORK_DIR}/qc/cutadapt
fastqc cutadapt/${baseName}.fastq.gz -t $nThreads -o ${WORK_DIR}/qc/cutadapt

		
#######################################################
## Align to genome with STAR                         ##
#######################################################


# align to genome
echo "aligning to genome..."
mkdir -p ${WORK_DIR}/bamSTAR
STAR --genomeDir ${GENOME_DIR}/sequence  --readFilesIn ${WORK_DIR}/cutadapt/${baseName}.fastq.gz --readFilesCommand zcat --outFileNamePrefix ${WORK_DIR}/bamSTAR/${baseName}_ --runThreadN $nThreads --outSAMmultNmax 1 --alignIntronMax 5000 --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --outWigType wiggle --outWigNorm RPM

samtools index ${WORK_DIR}/bamSTAR/${baseName}_Aligned.sortedByCoord.out.bam

${WIGtoBW_DIR}/wigToBigWig ${WORK_DIR}/bamSTAR/${baseName}_Signal.UniqueMultiple.str1.out.wig $chromSizesFile ${WORK_DIR}/bamSTAR/${baseName}_UniqueMultiple_F.bw
${WIGtoBW_DIR}/wigToBigWig ${WORK_DIR}/bamSTAR/${baseName}_Signal.UniqueMultiple.str2.out.wig $chromSizesFile ${WORK_DIR}/bamSTAR/${baseName}_UniqueMultiple_R.bw
${WIGtoBW_DIR}/wigToBigWig ${WORK_DIR}/bamSTAR/${baseName}_Signal.Unique.str1.out.wig $chromSizesFile ${WORK_DIR}/bamSTAR/${baseName}_Unique_F.bw
${WIGtoBW_DIR}/wigToBigWig ${WORK_DIR}/bamSTAR/${baseName}_Signal.Unique.str2.out.wig $chromSizesFile ${WORK_DIR}/bamSTAR/${baseName}_Unique_R.bw

rm ${WORK_DIR}/bamSTAR/${baseName}_Signal.UniqueMultiple.str1.out.wig
rm ${WORK_DIR}/bamSTAR/${baseName}_Signal.UniqueMultiple.str2.out.wig
rm ${WORK_DIR}/bamSTAR/${baseName}_Signal.Unique.str1.out.wig
rm ${WORK_DIR}/bamSTAR/${baseName}_Signal.Unique.str2.out.wig



## convert sam to bam, sort and index
#samtools view -@ $nThreads -b ${WORK_DIR}/bamSTAR/${baseName}_Aligned.out.sam | samtools sort -@ $nThreads -o ${WORK_DIR}/bamSTAR/${baseName}.sorted.bam -
#samtools index -@ $nThreads ${WORK_DIR}/bamSTAR/${baseName}.sorted.bam
#rm ${WORK_DIR}/bamSTAR/${baseName}_Aligned.out.sam


## create bigwig coverage tracks from bam
#mkdir -p ${WORK_DIR}/tracks
#Rscript makeCoverageBigwig.R ${WORK_DIR}/bamSTAR/${baseName}.sorted.bam ${WORK_DIR}/tracks/${baseName}_raw.bw raw
#Rscript makeCoverageBigwig.R ${WORK_DIR}/bamSTAR/${baseName}.sorted.bam ${WORK_DIR}/tracks/${baseName}_rpm.bw rpm


#######################################################
## Count reads with Salmon                           ##
#######################################################

# prepare reference genome
# See separate script indexTxpts4salmon.sh

######## NOTE: I do not have estimates for --fldMean and --fldSD as i have no access to the bioanalyser files #########
#######  therefore the quantification based on the effective transcript length will be wrong!!!! Not sure how important this ${WIGtoBW_DIR}#
# quantify mRNA transcripts
${SALMON_SING} salmon quant -i ${mRNAindex} -l A -r ${WORK_DIR}/cutadapt/${baseName}.fastq.gz --validateMappings -p ${nThreads} -o ${WORK_DIR}/salmon/mRNA/${baseName} --seqBias --gcBias --numBootstraps 100 




if [[ "${mRNAonly}"=="false" ]]
    then
    
    # quantify ncRNA transcripts
    ${SALMON_SING} salmon quant -i ${ncRNAindex} -l A -r ${WORK_DIR}/cutadapt/${baseName}.fastq.gz --validateMappings -p ${nThreads} -o ${WORK_DIR}/salmon/ncRNA/${baseName} --seqBias --gcBias --numBootstraps 100  
     
    # quantify pseudoRNA transcripts
    ${SALMON_SING} salmon quant -i ${pseudoIndex} -l A -r ${WORK_DIR}/cutadapt/${baseName}.fastq.gz --validateMappings -p ${nThreads} -o ${WORK_DIR}/salmon/pseudoRNA/${baseName} --seqBias --gcBias --numBootstraps 100  
     
    # quantify TnRNA transcripts
    ${SALMON_SING} salmon quant -i ${tnIndex} -l A -r ${WORK_DIR}/cutadapt/${baseName}.fastq.gz --validateMappings -p ${nThreads} -o ${WORK_DIR}/salmon/tnRNA/${baseName} --seqBias --gcBias --numBootstraps 100 
    
    
    kmerSize=15
    rptIndex=${GENOME_DIR}/sequence/${genomeVer}_repeats_index_${kmerSize}
    # quantify repeat transcripts
    ${SALMON_SING} salmon quant -i ${rptIndex} -l A -r ${WORK_DIR}/cutadapt/${baseName}.fastq.gz --validateMappings -p ${nThreads} -o ${WORK_DIR}/salmon/rptRNA_${kmerSize}/${baseName} --seqBias --gcBias --numBootstraps 100
    
    #######################################################
    ## Align to genome with STAR                         ##
    #######################################################
    
    
    # align to genome
    echo "aligning to genome..."
    mkdir -p ${WORK_DIR}/bamSTARrpts
    STAR --genomeDir ${GENOME_DIR}/sequence/repeats  --readFilesIn ${WORK_DIR}/cutadapt/${baseName}.fastq.gz --readFilesCommand zcat --outFileNamePrefix ${WORK_DIR}/bamSTARrpts/rpts_${baseName}_ --runThreadN $nThreads --outFilterMultimapNmax 6000 --outSAMmultNmax 1 --outFilterMismatchNmax 3 --alignIntronMax 1 --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --winAnchorMultimapNmax 6000 --outWigType wiggle --outWigNorm RPM 
    
    samtools index ${WORK_DIR}/bamSTARrpts/rpts_${baseName}_Aligned.sortedByCoord.out.bam

    ${WIGtoBW_DIR}/wigToBigWig ${WORK_DIR}/bamSTARrpts/rpts_${baseName}_Signal.UniqueMultiple.str1.out.wig $chromSizesFile ${WORK_DIR}/bamSTARrpts/rpts_${baseName}_UniqueMultiple_F.bw
    ${WIGtoBW_DIR}/wigToBigWig ${WORK_DIR}/bamSTARrpts/rpts_${baseName}_Signal.UniqueMultiple.str2.out.wig $chromSizesFile ${WORK_DIR}/bamSTARrpts/rpts_${baseName}_UniqueMultiple_R.bw
    ${WIGtoBW_DIR}/wigToBigWig ${WORK_DIR}/bamSTARrpts/rpts_${baseName}_Signal.Unique.str1.out.wig $chromSizesFile ${WORK_DIR}/bamSTARrpts/rpts_${baseName}_Unique_F.bw
    ${WIGtoBW_DIR}/wigToBigWig ${WORK_DIR}/bamSTARrpts/rpts_${baseName}_Signal.Unique.str2.out.wig $chromSizesFile ${WORK_DIR}/bamSTARrpts/rpts_${baseName}_Unique_R.bw
    
    rm ${WORK_DIR}/bamSTARrpts/rpts_${baseName}_Signal.UniqueMultiple.str1.out.wig
    rm ${WORK_DIR}/bamSTARrpts/rpts_${baseName}_Signal.UniqueMultiple.str2.out.wig
    rm ${WORK_DIR}/bamSTARrpts/rpts_${baseName}_Signal.Unique.str1.out.wig
    rm ${WORK_DIR}/bamSTARrpts/rpts_${baseName}_Signal.Unique.str2.out.wig
    
    # convert sam to bam, sort and index
    #samtools view -@ $nThreads -b ${WORK_DIR}/bamSTAR/rpts_${baseName}_Aligned.out.sam -o ${WORK_DIR}/bamSTAR/rpts_${baseName}.sorted.bam
    #| samtools sort -@ $nThreads -o ${WORK_DIR}/bamSTAR/rpts_${baseName}.sorted.bam -
    #samtools index -@ $nThreads ${WORK_DIR}/bamSTAR/rpts_${baseName}.sorted.bam
    #rm ${WORK_DIR}/bamSTAR/${baseName}_Aligned.out.sam
    
    
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
    annotFile_rpt=${GENOME_DIR}/annotation/c_elegans.PRJNA13758.${genomeVer}.annotations_rpt.gtf
    
    # count reads per feature with htseq
    mkdir -p ${WORK_DIR}/htseq
    htseq-count -f bam -r name -a 0 -m union --nonunique random  ${WORK_DIR}/bamBWA/${baseName}.bam $annotFile_rpt > ${WORK_DIR}/htseq/${baseName}_union_random.txt
    htseq-count -f bam -r name -a 0 -m union --nonunique none  ${WORK_DIR}/bamBWA/${baseName}.bam $annotFile_rpt > ${WORK_DIR}/htseq/${baseName}_union_none.txt
    #--additional-attr=gene_name
    
    samtools sort -T ${WORK_DIR}/bamBWA/${baseName}  -@ $nThreads -o ${WORK_DIR}/bamBWA/${baseName}_sort.bam  ${WORK_DIR}/bamBWA/${baseName}.bam
    samtools index -@ $nThreads ${WORK_DIR}/bamBWA/${baseName}_sort.bam
    rm ${WORK_DIR}/bamBWA/${baseName}.bam
fi

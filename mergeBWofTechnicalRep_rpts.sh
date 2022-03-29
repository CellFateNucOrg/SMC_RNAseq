#! /bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="mergebw"
#SBATCH --time=0-02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G

source $CONDA_ACTIVATE RNAseq

WORK_DIR=$PWD

genomeVer=WS275	
GENOME_DIR=${HOME}/genomeVer/${genomeVer}
chromSizesFile=$GENOME_DIR/annotation/ws235.chrom.sizes

echo "current date"
date

echo "current git branch"
git branch

echo "current git version"
git log -1 --format="%H"

mkdir -p ${WORK_DIR}/tracks_rpts

techRep=( `awk 'NR!=1 && NF>1 {print $3 "_" $4}' fastqList.txt | sort -u` )
biolRep=( `awk 'NR!=1 && NF>1 {print $3 "_" substr($4,1,length($4)-1)}' fastqList.txt | sort -u`)

# average technical reps for each strand
for i in ${techRep[@]}
do
        echo $i
        # technical reps forward strand
        repFiles=( `ls ${WORK_DIR}/bamSTARrpts/*_ribo0_UniqueMultiple_F.bw | grep ${i}` )
	echo merging technical reps F strand: ${repFiles[@]}
	wiggletools mean ${repFiles[@]} > ${WORK_DIR}/tracks_rpts/${i}_F_ribo0rpts_UniqueMultiple_RPM.wig
        wigToBigWig ${WORK_DIR}/tracks_rpts/${i}_F_ribo0rpts_UniqueMultiple_RPM.wig $chromSizesFile ${WORK_DIR}/tracks_rpts/${i}_F_ribo0rpts_UniqueMultiple_RPM.bw

        if [ -f ${WORK_DIR}/tracks_rpts/${i}_F_ribo0rpts_UniqueMultiple_RPM.bw ]; then
                rm ${WORK_DIR}/tracks_rpts/${i}_F_ribo0rpts_UniqueMultiple_RPM.wig
        fi

        # technical reps reverse strand
        repFiles=( `ls ${WORK_DIR}/bamSTARrpts/*_ribo0_UniqueMultiple_R.bw | grep ${i}` )
	echo merging technical reps R strand: ${repFiles[@]}
        wiggletools mean ${repFiles[@]} > ${WORK_DIR}/tracks_rpts/${i}_R_ribo0rpts_UniqueMultiple_RPM.wig
        wigToBigWig ${WORK_DIR}/tracks_rpts/${i}_R_ribo0rpts_UniqueMultiple_RPM.wig $chromSizesFile ${WORK_DIR}/tracks_rpts/${i}_R_ribo0rpts_UniqueMultiple_RPM.bw

        if [ -f ${WORK_DIR}/tracks_rpts/${i}_R_ribo0rpts_UniqueMultiple_RPM.bw ]; then
                rm ${WORK_DIR}/tracks_rpts/${i}_R_ribo0rpts_UniqueMultiple_RPM.wig
        fi
done


# average biological reps for each strand
for i in ${biolRep[@]}
do
	# biological reps forward strand
        repFiles=( `ls ${WORK_DIR}/tracks_rpts/*_ribo0rpts_UniqueMultiple_RPM.bw | grep ${i}` )
        echo merging biological reps F strand: ${repFiles[@]}
	wiggletools mean ${repFiles[@]} > ${WORK_DIR}/tracks_rpts/avr_${i}_F_ribo0rpts_UniqueMultiple_RPM.wig
        wigToBigWig ${WORK_DIR}/tracks_rpts/avr_${i}_F_ribo0rpts_UniqueMultiple_RPM.wig $chromSizesFile ${WORK_DIR}/tracks_rpts/avr_${i}_F_ribo0rpts_UniqueMultiple_RPM.bw

        if [ -f ${WORK_DIR}/tracks_rpts/avr_${i}_F_ribo0rpts_UniqueMultiple_RPM.bw ]; then
                rm ${WORK_DIR}/tracks_rpts/avr_${i}_F_ribo0rpts_UniqueMultiple_RPM.wig
        fi

	# biological reps reverse strand
        repFiles=( `ls ${WORK_DIR}/tracks_rpts/*_ribo0rpts_UniqueMultiple_RPM.bw | grep ${i}` )
	echo merging biological reps R strand: ${repFiles[@]}
        wiggletools mean ${repFiles[@]} > ${WORK_DIR}/tracks_rpts/avr_${i}_R_ribo0rpts_UniqueMultiple_RPM.wig
        wigToBigWig ${WORK_DIR}/tracks_rpts/avr_${i}_R_ribo0rpts_UniqueMultiple_RPM.wig $chromSizesFile ${WORK_DIR}/tracks_rpts/avr_${i}_R_ribo0rpts_UniqueMultiple_RPM.bw

        if [ -f ${WORK_DIR}/tracks_rpts/avr_${i}_R_ribo0rpts_UniqueMultiple_RPM.bw ]; then
                rm ${WORK_DIR}/tracks_rpts/avr_${i}_R_ribo0rpts_UniqueMultiple_RPM.wig
        fi
done

# sum forward and reverse strand
for i in ${biolRep[@]}
do
	repFiles=( `ls ${WORK_DIR}/tracks_rpts/avr_*_ribo0rpts_UniqueMultiple_RPM.bw | grep ${i}` )
	echo summing F and R  strands: ${repFiles[@]}

	bigWigMerge ${repFiles[@]} ${WORK_DIR}/tracks_rpts/sumFR_${i}_ribo0rpts_UniqueMultiple.bedGraph

	bedGraphToBigWig ${WORK_DIR}/tracks_rpts/sumFR_${i}_ribo0rpts_UniqueMultiple.bedGraph $chromSizesFile ${WORK_DIR}/tracks_rpts/sumFR_${i}_ribo0rpts_UniqueMultiple.bw

	if [ -f ${WORK_DIR}/tracks_rpts/sumFR_${i}_ribo0rpts_UniqueMultiple.bw ]; then
		rm ${WORK_DIR}/tracks_rpts/sumFR_${i}_ribo0rpts_UniqueMultiple.bedGraph
	fi
done

########### make ratios
Rscript makeBWratios_rpts.R
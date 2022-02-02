#! /bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="RNAseq"
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

mkdir -p ${WORK_DIR}/tracks

biolRep=( `awk 'NR!=1 && NF>1 {print $3 "_" $4}' fastqList.txt | sort -u` )

for i in ${biolRep[@]}
do
	echo $i
	repFiles=( `ls ${WORK_DIR}/bamSTAR/*_UniqueMultiple.bw | grep $i` )

	bigWigMerge ${repFiles[@]} ${WORK_DIR}/tracks/${i}_STAR_UniqueMultiple.bedGraph

	bedGraphToBigWig ${WORK_DIR}/tracks/${i}_STAR_UniqueMultiple.bedGraph $chromSizesFile ${WORK_DIR}/tracks/${i}_STAR_UniqueMultiple.bw

	if [ -f ${WORK_DIR}/tracks/STAR_${i}.bw ]; then
		rm ${WORK_DIR}/tracks/${i}_STAR_UniqueMultiple.bedGraph
	fi
done




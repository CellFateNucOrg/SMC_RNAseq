#! /bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="RNAseq"
##SBATCH --partition=epyc2
##SBATCH --qos=job_epyc2
#SBATCH --time=0-02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1#2-9,11-83%12

source $CONDA_ACTIVATE RNAseq

echo "current date"
date

echo "current git branch"
git branch

echo "current git version"
git log -1 --format="%H"

# check if paired end or single end
grep "\sfileName2\s" fastqList.txt
if [ "$?" == 0  ]
then
        isPE=true
else
        echo "running single end commands"
fi

fastqFileList=./fastqList.txt
if [ "$isPE" == "true" ]; then
        fastqFile1s=(`cut -f1 $fastqFileList`)
        fastqFile2s=(`cut -f2 $fastqFileList`)
        sampleNames=(`cut -f3 $fastqFileList`)
        repeatNums=(`cut -f4 $fastqFileList`)
        laneNums=(`cut -f5 $fastqFileList`)
else
        fastqFile1s=(`cut -f1 $fastqFileList`)
        sampleNames=(`cut -f2 $fastqFileList`)
        repeatNums=(`cut -f3 $fastqFileList`)
        laneNums=(`cut -f4 $fastqFileList`)
fi



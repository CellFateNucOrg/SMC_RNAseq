#! /bin/bash

DIRNAME=$1

if [ -z $DIRNAME ]; then
  echo "Provide name of fastqc results directory to process"
  exit
fi

fqcFiles=( `ls $DIRNAME/*.zip` )

if [[ -f "${DIRNAME}/readCount.txt" ]]; then
  rm ${DIRNAME}/readCount.txt
fi
touch ${DIRNAME}/readCount.txt
#f=${fqcFiles[0]}
for f in ${fqcFiles[@]}
do
	unzip -d $DIRNAME $f
	fqcDIR=`basename $f`
	fqcDIR=${fqcDIR%.zip}
	qcDIR=`basename $DIRNAME`
	totalSeq=`grep "Total Sequences" ${DIRNAME}/${fqcDIR}/fastqc_data.txt`
	totalSeq=`echo ${totalSeq#Total Sequences} | xargs` # use xargs to remove leading and trailing spaces
	echo ${qcDIR}$'\t'${fqcDIR}$'\t'${totalSeq#Total Sequences} >> ${DIRNAME}/readCount.txt
	rm -r $DIRNAME/$fqcDIR/*
	rmdir $DIRNAME/$fqcDIR
done

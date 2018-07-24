#!/bin/bash
BOWTIE=$1
SAMTOOLS=$2
INFILE=$3
OUTFILE=$4
REFFILE=$5
TEMPPREFIX=$6
BOWTIEARGS="${@:7}"

FQFILE=$INFILE
if [ ${INFILE: -4} == ".bam" ]
then
	echo "Step: Convert bam to fastq"
	FQFILE=${TEMPPREFIX}_tmp01.fq
	bedtools bamtofastq -i ${INFILE} -fq ${FQFILE}
fi

echo "Step: Alignment"
${BOWTIE} ${BOWTIEARGS} -x ${REFFILE} -U ${FQFILE} -S ${TEMPPREFIX}_tmp02.sam
echo "Step: Converting to bam"
${SAMTOOLS} view -b -T ${REFFILE} -o ${TEMPPREFIX}_tmp03_unsorted.bam ${TEMPPREFIX}_tmp02.sam
echo "Step: Sorting bam"
${SAMTOOLS} sort -m 4G -o ${OUTFILE} ${TEMPPREFIX}_tmp03_unsorted.bam
echo "Step: Indexing bam"
${SAMTOOLS} index ${OUTFILE}

# rm ${TEMPPREFIX}_tmp01.sai
# rm ${TEMPPREFIX}_tmp02.sam
# rm ${TEMPPREFIX}_tmp03_unsorted.bam

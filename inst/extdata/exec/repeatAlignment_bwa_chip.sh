#!/bin/bash
BWA=$1
SAMTOOLS=$2
INFILE=$3
OUTFILE=$4
REFFILE=$5
TEMPPREFIX=$6
BWAARGS="${@:7}"

echo "Step: aln"
if [ ${INFILE: -4} == ".bam" ]
then
	${BWA} aln ${BWAARGS} -b ${REFFILE} ${INFILE} > ${TEMPPREFIX}_tmp01.sai
else
	${BWA} aln ${BWAARGS} ${REFFILE} ${INFILE} > ${TEMPPREFIX}_tmp01.sai
fi

echo "Step: samse"
${BWA} samse ${REFFILE} ${TEMPPREFIX}_tmp01.sai ${INFILE} > ${TEMPPREFIX}_tmp02.sam
echo "Step: Converting to bam"
${SAMTOOLS} view -b -T ${REFFILE} -o ${TEMPPREFIX}_tmp03_unsorted.bam ${TEMPPREFIX}_tmp02.sam
echo "Step: Sorting bam"
${SAMTOOLS} sort -m 4G -o ${OUTFILE} ${TEMPPREFIX}_tmp03_unsorted.bam
echo "Step: Indexing bam"
${SAMTOOLS} index ${OUTFILE}

# rm ${TEMPPREFIX}_tmp01.sai
# rm ${TEMPPREFIX}_tmp02.sam
# rm ${TEMPPREFIX}_tmp03_unsorted.bam

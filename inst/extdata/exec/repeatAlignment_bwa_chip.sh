#!/bin/bash
BWA=$1
SAMTOOLS=$2
INFILE=$3
OUTFILE=$4
REFFILE=$5
TEMPPREFIX=$6
BWAARGS="${@:7}"

${BWA} aln ${BWAARGS} ${REFFILE} ${INFILE} > ${TEMPPREFIX}_tmp01.sai
${BWA} samse ${REFFILE} ${TEMPPREFIX}_tmp01.sai ${INFILE} > ${TEMPPREFIX}_tmp02.sam
${SAMTOOLS} view -b -T ${REFFILE} -o ${TEMPPREFIX}_tmp03_unsorted.bam ${TEMPPREFIX}_tmp02.sam
${SAMTOOLS} sort -o ${OUTFILE} ${TEMPPREFIX}_tmp03_unsorted.bam
${SAMTOOLS} index ${OUTFILE}

# rm ${TEMPPREFIX}_tmp01.sai
# rm ${TEMPPREFIX}_tmp02.sam
# rm ${TEMPPREFIX}_tmp03_unsorted.bam

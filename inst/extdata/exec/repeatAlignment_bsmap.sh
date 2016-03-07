#!/bin/bash
BSMAP=$1
SAMTOOLS=$2
INFILE=$3
OUTFILE=$4
REFFILE=$5
TEMPPREFIX=$6
BSMAPARGS="${@:7}"

# echo "${BSMAP} -a ${INFILE} -d ${REFFILE} -o ${TEMPPREFIX}_unsorted.bam ${BSMAPARGS}"
${BSMAP} -a ${INFILE} -d ${REFFILE} -o ${TEMPPREFIX}_unsorted.bam ${BSMAPARGS}
${SAMTOOLS} sort -o ${OUTFILE} ${TEMPPREFIX}_unsorted.bam
${SAMTOOLS} index ${OUTFILE}

# rm ${TEMPPREFIX}_unsorted.bam

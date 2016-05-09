#!/bin/bash
BOWTIE=$1
SAMTOOLS=$2
INFILE=$3
OUTFILE=$4
REFFILE=$5
TEMPPREFIX=$6
BOWTIEARGS="${@:7}"

echo "Step: Convert bam to fastq"
bedtools bamtofastq -i ${INFILE} -fq ${TEMPPREFIX}_tmp01.fq
echo "Step: Alignment"
${BOWTIE} ${BOWTIEARGS} -x ${REFFILE} -U ${TEMPPREFIX}_tmp01.fq -S ${TEMPPREFIX}_tmp02.sam
echo "Step: Converting to bam"
${SAMTOOLS} view -b -T ${REFFILE} -o ${TEMPPREFIX}_tmp03_unsorted.bam ${TEMPPREFIX}_tmp02.sam
echo "Step: Sorting bam"
${SAMTOOLS} sort -m 4G -o ${OUTFILE} ${TEMPPREFIX}_tmp03_unsorted.bam
echo "Step: Indexing bam"
${SAMTOOLS} index ${OUTFILE}

# rm ${TEMPPREFIX}_tmp01.sai
# rm ${TEMPPREFIX}_tmp02.sam
# rm ${TEMPPREFIX}_tmp03_unsorted.bam

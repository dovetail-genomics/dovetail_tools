#!/usr/bin/env bash

ref=$1
bam=$2
peaks=$3
prefix=$4

samtools faidx ${ref}
cut -f1,2 ${ref}".fai" > ${prefix}".genome"
genome=${prefix}".genome"

SRCDIR=`dirname $0`

#rerorder peaks file based on order of chromosomes in reference

echo "bedtools sort -g ${genome} -i ${peaks} > ${prefix}_reordered_peaks.bed"
bedtools sort -g ${genome} -i ${peaks} > ${prefix}_reordered_peaks.bed

bed=${prefix}"_reordered_peaks.bed"

#compute how many reads intersect with peaks
bedtools intersect -a ${bam} -b ${bed} -bed | sort -k4 > ${prefix}_peak_intersect.bed &
bedtools window -w 500 -abam ${bam} -b ${peaks} -bed | sort -k4 > ${prefix}_peaks_intersect_500.bed &
bedtools window -w 1000 -abam ${bam} -b ${peaks} -bed | sort -k4 > ${prefix}_peaks_intersect_1000.bed &
bedtools window -w 2000 -abam ${bam} -b ${peaks} -bed | sort -k4 > ${prefix}_peaks_intersect_2000.bed &
bedtools window -w 5000 -abam ${bam} -b ${peaks} -bed | sort -k4 > ${prefix}_peaks_intersect_5000.bed &

wait

#generate peak enrichment plot

bamCoverage --bam ${bam} --outFileName  ${prefix}_coverage.bigwig --outFileFormat bigwig
computeMatrix reference-point --scoreFileName ${prefix}_coverage.bigwig --referencePoint \
 	center --beforeRegionStartLength 10000 --afterRegionStartLength 10000 \
 	--outFileName ${prefix}_coverage.matrix.gz

plotProfile -m ${prefix}_coverage.matrix.gz --perGroup  -out ${prefix}_peak_enrichment.png


#print final stats
python ${SRCDIR}/count.py -b1  ${prefix}_peak_intersect.bed -b2 ${prefix}_peaks_intersect_500.bed \
	-b3 ${prefix}_peaks_intersect_1000.bed -b4 ${prefix}_peaks_intersect_2000.bed \
	-b5 ${prefix}_peaks_intersect_5000.bed   -bam ${bam} -peaks ${peaks}



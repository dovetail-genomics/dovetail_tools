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

bedtools sort -g ${genome} -i ${bed} > ${prefix}_reordered_peaks.bed

bed=${prefix}"_reordered_peaks.bed"

#compute how many reads intersect with peaks
bedtools intersect -a ${bam} -b {bed} | sort -k4 > ${prefix}_peak_intersect.bed &
bedtools window -w 500 -abam ${bam} -b ${peaks} -bed | sort -k4 > ${prefix}_peaks_intersect_500.bed &
bedtools window -w 1000 -abam ${bam} -b ${peaks} -bed | sort -k4 > ${prefix}_peaks_intersect_1000.bed &
bedtools window -w 2000 -abam ${bam} -b ${peaks} -bed | sort -k4 > ${prefix}_peaks_intersect_2000.bed &
bedtools window -w 5000 -abam ${bam} -b ${peaks} -bed | sort -k4 > ${prefix}_peaks_intersect_5000.bed &

wait


#widen peaks file to accomodate reads around peaks
bedtools slop -b 1000 -i ${peaks} -g ${genome} > ${prefix}_peaks_1000_wide.bed &
bedtools slop -b 2000 -i ${peaks} -g ${genome} > ${prefix}_peaks_2000_wide.bed &
bedtools slop -b 5000 -i ${peaks} -g ${genome} > ${prefix}_peaks_5000_wide.bed &

wait

#complement peaks for other regions
bedtools complement -i ${bed} -g ${genome} > ${prefix}_peak_complement.bed &
bedtools complement -i ${prefix}_peaks_1000_wide.bed -g ${genome} > ${prefix}_peaks_1000_wide_complement.bed &
bedtools complement -i ${prefix}_peaks_2000_wide.bed -g ${genome} > ${prefix}_peaks_2000_wide_complement.bed &
bedtools complement -i ${prefix}_peaks_5000_wide.bed -g ${genome} > ${prefix}_peaks_5000_wide_complement.bed &

wait



#compute read coverage in and outside peaks and widened regions

bedtools coverage -sorted -g {genome} -a ${bed} -b ${bam} -mean > ${prefix}.peak.coverage.bedgraph &
bedtools coverage -sorted -g ${genome} -a ${prefix}_peak_complement.bed -b ${bam} -mean > ${prefix}.other.coverage.bedgraph &

wait

bedtools coverage -sorted -g {genome} -a $${prefix}_peaks_1000_wide.bed-b ${bam} -mean > ${prefix}.peak.coverage.1000.bedgraph &
bedtools coverage -sorted -g ${genome} -a {prefix}_peaks_1000_wide_complement.bed -b ${bam} -mean > ${prefix}.other.coverage.1000.bedgraph &
wait

bedtools coverage -sorted -g {genome} -a $${prefix}_peaks_2000_wide.bed-b ${bam} -mean > ${prefix}.peak.coverage.2000.bedgraph &
bedtools coverage -sorted -g ${genome} -a {prefix}_peaks_2000_wide_complement.bed -b ${bam} -mean > ${prefix}.other.coverage.2000.bedgraph &
wait

bedtools coverage -sorted -g {genome} -a $${prefix}_peaks_5000_wide.bed-b ${bam} -mean > ${prefix}.peak.coverage.5000.bedgraph &
bedtools coverage -sorted -g ${genome} -a {prefix}_peaks_5000_wide_complement.bed -b ${bam} -mean > ${prefix}.other.coverage.5000.bedgraph &
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
	-b5 ${prefix}_peaks_intersect_5000.bed -cp1 ${prefix}.peak.coverage.bedgraph \
	-co1 {prefix}.other.coverage.bedgraph -cp2 ${prefix}.peak.coverage.1000.bedgraph \
	-co2 ${prefix}.other.coverage.1000.bedgraph -cp3 ${prefix}.peak.coverage.2000.bedgraph \
	-co3 ${prefix}.other.coverage.2000.bedgraph -cp4 ${prefix}.peak.coverage.5000.bedgraph \
	-co4 ${prefix}.other.coverage.5000.bedgraph



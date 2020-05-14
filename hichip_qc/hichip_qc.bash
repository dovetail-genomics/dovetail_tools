#!/usr/bin/env bash


if [  $# -le 3 ]
then
        echo "Too few arguments. Please provide all the required arguments."
        echo "Usage: ./hichip_qc.bash <reference_fasta> <alignment_bam> <chipseq_peaks>  <output_prefix>"
        exit 1
fi

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

bed_chr20=${prefix}"_chr20_reordered_peaks.bed"
grep -w 'chr20' ${bed} > ${bed_chr20}

#compute how many reads intersect with peaks
bedtools intersect -a ${bam} -b ${bed} -bed | sort -k4 > ${prefix}_peak_intersect.bed &
bedtools window -w 500 -abam ${bam} -b ${peaks} -bed | sort -k4 > ${prefix}_peaks_intersect_500.bed &
bedtools window -w 1000 -abam ${bam} -b ${peaks} -bed | sort -k4 > ${prefix}_peaks_intersect_1000.bed &
bedtools window -w 2000 -abam ${bam} -b ${peaks} -bed | sort -k4 > ${prefix}_peaks_intersect_2000.bed &
bedtools window -w 5000 -abam ${bam} -b ${peaks} -bed | sort -k4 > ${prefix}_peaks_intersect_5000.bed &

wait

#generate peak enrichment plot

bamCoverage --bam ${bam} --outFileName  ${prefix}_coverage.bigwig --outFileFormat bigwig & 

plotFingerprint -b ${bam} --region chr20 --plotFile ${prefix}_chip_fingerprint_plot.png --outRawCounts ${prefix}_counts.tab &

wait 

python ${SRCDIR}/plot_chip_fingerprint.py -table ${prefix}_counts.tab -output ${prefix}_chip_fingerprint_plot.png 
python ${SRCDIR}/plot_chip_enrichment.py -bam ${bam} -peaks ${peaks} -output ${prefix}_chip_enrichment_plot.png

#print final stats
python ${SRCDIR}/count.py -b1  ${prefix}_peak_intersect.bed -b2 ${prefix}_peaks_intersect_500.bed \
	-b3 ${prefix}_peaks_intersect_1000.bed -b4 ${prefix}_peaks_intersect_2000.bed \
	-b5 ${prefix}_peaks_intersect_5000.bed   -bam ${bam} -peaks ${peaks}



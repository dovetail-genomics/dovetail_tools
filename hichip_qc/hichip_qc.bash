#!/usr/bin/env bash


if [  $# -le 5 ]
then
        echo "Too few arguments. Please provide all the required arguments."
        echo "Usage: ./hichip_qc.bash <reference_fasta> <read1_fastq> <read2_fastq>  <chipseq_peaks>  <output_prefix> <num_cores>"
        exit 1
fi

ref=$1
r1fq=$2
r2fq=$3
peaks=$4
prefix=$5
cores=$6

sample=`basename ${prefix}`
SRCDIR=`dirname $0`

OUTPUTFILE=${prefix}"_hichip_qc_metrics.txt"
TMPOUT=${prefix}"_hichip_qc_metrics.txt.tmp"
#first run omnic_qc

bam=${prefix}"-PT.bam"


genome=${prefix}.genome
${SRCDIR}/../omni-c_qc.bash ${ref} ${r1fq} ${r2fq} ${prefix}  ${sample} ${cores} > ${TMPOUT}



bedtools sort -g ${genome} -i ${peaks} > ${prefix}_reordered_peaks.bed

bed=${prefix}"_reordered_peaks.bed"

bed_chr20=${prefix}"_chr20_reordered_peaks.bed"
grep -w 'chr20' ${bed} > ${bed_chr20}

#find reads in the blacklisted region

bedtools intersect -a ${bam} -b ${SRCDIR}/hg38.blacklist.bed -bed | sort -k4 > ${prefix}_blacklist_intersect.bed 


#compute how many reads intersect with peaks
awk '{$2=$2+$10-50;print}' ${peaks} | awk '{$3=$2+100;print}' | awk '{print $1"\t"$2"\t"$3}' > ${prefix}_tmp.bed
bedtools intersect -a ${bam} -b ${prefix}_tmp.bed -bed | sort -k4 > ${prefix}_peak_intersect_100.bed
awk '{$2=$2+$10-100;print}' ${peaks} | awk '{$3=$2+200;print}' |  awk '{print $1"\t"$2"\t"$3}'> ${prefix}_tmp.bed
bedtools intersect -a ${bam} -b ${prefix}_tmp.bed -bed | sort -k4 | sort -k4 > ${prefix}_peaks_intersect_200.bed 
awk '{$2=$2+$10-250;print}' ${peaks} | awk '{$3=$2+500;print}' |  awk '{print $1"\t"$2"\t"$3}' > ${prefix}_tmp.bed
bedtools intersect -a ${bam} -b ${prefix}_tmp.bed -bed | sort -k4 > ${prefix}_peaks_intersect_500.bed 


wait

#generate peak enrichment plot

bamCoverage --bam ${bam} --outFileName  ${prefix}_coverage.bigwig --outFileFormat bigwig -p ${cores}& 

#plotFingerprint -b ${bam} --region chr20 --plotFile ${prefix}_chip_fingerprint_plot.png --outRawCounts ${prefix}_counts.tab &


wait 

python ${SRCDIR}/plot_chip_enrichment.py -bam ${bam} -peaks ${peaks} -output ${prefix}_chip_enrichment_plot.png

python ${SRCDIR}/count.py -b1 ${prefix}_peak_intersect_100.bed  -b2 ${prefix}_peaks_intersect_200.bed \
	-b3  ${prefix}_peaks_intersect_500.bed  -b4 ${prefix}_blacklist_intersect.bed -bam ${bam} -peaks ${peaks} >> ${TMPOUT}

mv  ${TMPOUT} $OUTPUTFILE


rm ${prefix}_reordered_peaks.bed  ${bed_chr20} ${prefix}_peak_intersect_100.bed ${prefix}_peaks_intersect_200.bed  \
	${prefix}_peaks_intersect_500.bed ${TMPOUT} ${prefix}_tmp.bed

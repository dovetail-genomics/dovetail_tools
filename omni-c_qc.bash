#!/usr/bin/env bash


if [  $# -le 5 ]
then
        echo "Too few arguments. Please provide all the required arguments."
        echo "Usage: ./omni-c_qc.bash <reference_fasta> <read1_fastq> <reaf2_fastq>  <output_prefix>  <sample_name> <num_cores>"
        exit 1
fi




ref=$1
fq1=$2
fq2=$3
outprefix=$4
rg=$5
cores=$6

#get source directory
SRCDIR=`dirname $0`

samtools faidx ${ref}
cut -f1,2 ${ref}".fai" > ${outprefix}".genome"
genome=${outprefix}".genome"

bwa mem -5SP -T0 -t${cores} \
    -R "@RG\tID:$rg\tSM:$rg\tLB:$rg\tPL:ILLUMINA\tPU:none" \
    $ref \
    $fq1 \
    $fq2 \
| pairtools parse \
	--chroms-path ${genome} \
	--min-mapq 40 \
	--walks-policy 5unique \
	--max-inter-align-gap 30 \
      	--nproc-in ${cores} --nproc-out ${cores} \
| pairtools sort  \
	--nproc ${cores} \
| pairtools dedup \
	--nproc-in ${cores} --nproc-out ${cores} \
	--mark-dups  \
 	--output-stats ${outprefix}-PT.stats.txt \
	--output-dups - \
| pairtools split \
	--nproc-in ${cores} --nproc-out ${cores} \
	--output-pairs ${outprefix}.PT.pairs.gz  \
	--output-sam - \
| samtools view -bS - \
| samtools sort -@${cores} - -o ${outprefix}-PT.bam


samtools index ${outprefix}-PT.bam;


preseq lc_extrap -B -P -e 2.1e9 -s 1e8 -seg_len 1000000000 -o $outprefix.preseq ${outprefix}-PT.bam

${SRCDIR}/get_qc.py -p ${outprefix}-PT.stats.txt -d ${outprefix}.preseq


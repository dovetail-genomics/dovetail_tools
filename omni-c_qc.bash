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
| pairtools parse --chroms-path ${genome}\
| pairtools sort  --nproc ${cores} \
| pairtools dedup --nproc-in ${cores} --nproc-out ${cores} --mark-dups  \
 	--output-stats ${outprefix}-PT.stats.txt --output-dups - \
| pairtools split --nproc-in ${cores} --nproc-out ${cores} \
	--output-pairs ${outprefix}.PT.pairs.gz  \
	--output-sam - \
| ${SRCDIR}/add_mate_MQ.py \
| samtools view -bS - \
| samtools sort -@${cores} - -o ${outprefix}-PT.bam


samtools index ${outprefix}-PT.bam;


preseq lc_extrap -B -P -e 2.1e9 -s 1e8 -seg_len 1000000000 -o $outprefix.preseq ${outprefix}-PT.bam

ps300m=`cat ${outprefix}.preseq | grep -P "^300000000.0" | awk '{print $2}'`

qualthresh=40
mate_filter_cmd="perl -e 'while (<STDIN>) { m/MQ:i:(\\d+)/; if (\$1 >= \$ARGV[0]) { print; }}' $qualthresh"

r1=`samtools view -c -q $qualthresh -f 0x40 -F 2304 ${outprefix}-PT.bam`
r2=`samtools view -c -q $qualthresh -f 0x80 -F 2304 ${outprefix}-PT.bam`

mapped_pairs=`samtools view -q $qualthresh -f 0x40 -F 2316 ${outprefix}-PT.bam | eval $mate_filter_cmd | wc -l`
pcr_dupe_pairs=`samtools view -q $qualthresh -u -f 0x40 -F 2316 ${outprefix}-PT.bam | samtools view -f 0x400 | eval $mate_filter_cmd | wc -l`

mapped_nondupe_pairs=`samtools view -q $qualthresh -f 0x40 -F 3340 ${outprefix}-PT.bam | eval $mate_filter_cmd | wc -l`
mapped_nondupe_pairs_cis=`samtools view -q $qualthresh -f 0x40 -F 3340 ${outprefix}-PT.bam | awk '{if (sqrt($9^2) > 0) { print; }}' | eval $mate_filter_cmd | wc -l`
mapped_nondupe_pairs_cis_lt1000=`samtools view -q $qualthresh -f 0x40 -F 3340 ${outprefix}-PT.bam | awk '{if ((sqrt($9^2) <= 1000) && ($9 != 0)) { print; }}' | eval $mate_filter_cmd | wc -l`
mapped_nondupe_pairs_cis_gt1000=`samtools view -q $qualthresh -f 0x40 -F 3340 ${outprefix}-PT.bam | awk '{if (sqrt($9^2) > 1000) { print; }}' | eval $mate_filter_cmd | wc -l`
mapped_nondupe_pairs_cis_gt10000=`samtools view -q $qualthresh -f 0x40 -F 3340 ${outprefix}-PT.bam | awk '{if (sqrt($9^2) > 10000) { print; }}' | eval $mate_filter_cmd | wc -l`
mapped_trans_pairs=`samtools view -q $qualthresh -f 0x40 -F 3340 ${outprefix}-PT.bam | awk '{if (sqrt($9^2) == 0) { print; }}' | eval $mate_filter_cmd | wc -l`

valid_pairs=$(($mapped_nondupe_pairs_cis_gt1000 + $mapped_trans_pairs))


echo "Mapping Quality Threshold         :" $qualthresh
echo "Read1                             :" `echo $r1 | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
echo "Read2                             :" `echo $r2 | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
echo "Mapped pairs                      :" `echo $mapped_pairs | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
echo "PCR dupe pairs                    :" `echo $pcr_dupe_pairs | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
echo "Mapped nondupe pairs              :" `echo $mapped_nondupe_pairs | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
echo "Valid Pairs (cis>1000bp + trans)  :" `echo $valid_pairs | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
echo "Mapped nondupe pairs cis          :" `echo $mapped_nondupe_pairs_cis | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
echo "Mapped nondupe pairs cis <=1000bp :" `echo $mapped_nondupe_pairs_cis_lt1000 | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
echo "Mapped nondupe pairs cis >1000bp  :" `echo $mapped_nondupe_pairs_cis_gt1000 | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
echo "Mapped nondupe pairs cis >10000bp :" `echo $mapped_nondupe_pairs_cis_gt10000 | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
echo "Mapped nondupe trans pairs        :" `echo $mapped_trans_pairs | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
echo "Expected unique pairs at 300M sequencing: " $ps300m

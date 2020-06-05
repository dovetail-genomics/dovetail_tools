#!/usr/bin/env bash

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
 	--output-stats ${outprefix}-PT.stats.txt \
| pairtools split --nproc-in ${cores} --nproc-out ${cores} \
	--output-pairs ${outprefix}.PT.pairs.gz  \
	--output-sam - \
| samtools view -bS - \
| samtools sort -@${cores} - -o ${outprefix}-PT.bam


samtools index ${outprefix}-PT.bam;


preseq lc_extrap -B -P -e 2.1e9 -s 1e8 -seg_len 1000000000 -o $outprefix.preseq ${outprefix}-PT.bam

ps300m=`cat ${outprefix}.preseq | grep -P "^300000000.0" | awk '{print $2}'`

qualthresh=40
mate_filter_cmd="perl -e 'while (<STDIN>) { m/MQ:i:(\\d+)/; if (\$1 >= \$ARGV[0]) { print; }}' $qualthresh"

r1=`samtools view -c -q $qualthresh -f 0x40 -F 2304 ${outprefix}-PT.bam`
r2=`samtools view -c -q $qualthresh -f 0x80 -F 2304 ${outprefix}-PT.bam`

mapped_pairs=`samtools view -q $qualthresh -f 0x40 -F 2316 ${outprefix}-PT.bam | perl -e 'while (<STDIN>) { m/MQ:i:(\\d+)/; if (\$1 >= \$ARGV[0]) { print; }}' | wc -l`
pcr_dupe_pairs=`samtools view -q $qualthresh -u -f 0x40 -F 2316 ${outprefix}-PT.bam | samtools view -f 0x400 | perl -e 'while (<STDIN>) { m/MQ:i:(\\d+)/; if (\$1 >= \$ARGV[0]) { print; }}' | wc -l`

mapped_nondupe_pairs=`samtools view -q $qualthresh -f 0x40 -F 3340 ${outprefix}-PT.bam | perl -e 'while (<STDIN>) { m/MQ:i:(\\d+)/; if (\$1 >= \$ARGV[0]) { print; }}' | wc -l`
mapped_nondupe_pairs_cis=`samtools view -q $qualthresh -f 0x40 -F 3340 ${outprefix}-PT.bam | awk '{if (sqrt($9^2) > 0) { print; }}' | perl -e 'while (<STDIN>) { m/MQ:i:(\\d+)/; if (\$1 >= \$ARGV[0]) { print; }}' | wc -l`
mapped_nondupe_pairs_cis_lt1000=`samtools view -q $qualthresh -f 0x40 -F 3340 ${outprefix}-PT.bam | awk '{if ((sqrt($9^2) <= 1000) && ($9 != 0)) { print; }}' | perl -e 'while (<STDIN>) { m/MQ:i:(\\d+)/; if (\$1 >= \$ARGV[0]) { print; }}' | wc -l`
mapped_nondupe_pairs_cis_gt1000=`samtools view -q $qualthresh -f 0x40 -F 3340 ${outprefix}-PT.bam | awk '{if (sqrt($9^2) > 1000) { print; }}' | perl -e 'while (<STDIN>) { m/MQ:i:(\\d+)/; if (\$1 >= \$ARGV[0]) { print; }}' | wc -l`
mapped_nondupe_pairs_cis_gt10000=`samtools view -q $qualthresh -f 0x40 -F 3340 ${outprefix}-PT.bam | awk '{if (sqrt($9^2) > 10000) { print; }}' | perl -e 'while (<STDIN>) { m/MQ:i:(\\d+)/; if (\$1 >= \$ARGV[0]) { print; }}' | wc -l`
mapped_trans_pairs=`samtools view -q $qualthresh -f 0x40 -F 3340 ${outprefix}-PT.bam | awk '{if (sqrt($9^2) == 0) { print; }}' | perl -e 'while (<STDIN>) { m/MQ:i:(\\d+)/; if (\$1 >= \$ARGV[0]) { print; }}' | wc -l`

valid_pairs=$(($mapped_nondupe_pairs_cis_gt1000 + $mapped_trans_pairs))

echo "Mapping Quality Threshold         :" $qualthresh
echo "Read1                             :" $r1
echo "Read2                             :" $r2
echo "Mapped pairs                      :" $mapped_pairs
echo "PCR dupe pairs                    :" $pcr_dupe_pairs
echo "Mapped nondupe pairs              :" $mapped_nondupe_pairs
echo "Valid Pairs (cis>1000bp + trans)  :" $valid_pairs
echo "Mapped nondupe pairs cis          :" $mapped_nondupe_pairs_cis
echo "Mapped nondupe pairs cis <=1000bp :" $mapped_nondupe_pairs_cis_lt1000
echo "Mapped nondupe pairs cis >1000bp  :" $mapped_nondupe_pairs_cis_gt1000
echo "Mapped nondupe pairs cis >10000bp :" $mapped_nondupe_pairs_cis_gt10000
echo "Mapped nondupe trans pairs        :" $mapped_trans_pairs
echo "Expected unique pairs at 300M sequencing: " $ps300m

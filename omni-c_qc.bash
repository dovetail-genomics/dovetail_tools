#!/usr/bin/env bash

ref=$1
fq1=$2
fq2=$3
out=$4
rg=$5

#get source directory
SRCDIR=`dirname $0`

bwa mem -5SP -T0 -t4 \
    -R "@RG\tID:$rg\tSM:$rg\tLB:$rg\tPL:ILLUMINA\tPU:none" \
    $ref \
    $fq1 \
    $fq2 \
| samblaster -i stdin -o stdout \
| $SRCDIR/add_mate_MQ.py \
| samtools view -S -h -b \
| samtools sort --threads 4 - > $out ;

samtools index $out ;

preseq lc_extrap -B -P -e 2.1e9 -s 1e8 -seg_len 1000000000 -o $out.preseq $out

ps300m=`cat $out.preseq | grep -P "^300000000.0" | awk '{print $2}'`

qualthresh=40
mate_filter_cmd="perl -ne 'm/MQ:i:(\\d+)/; if (\$1 >= \$ENV{qualthresh}) { print; }'"

r1=`samtools view -c -q $qualthresh -f 0x40 -F 2304 $out`
r2=`samtools view -c -q $qualthresh -f 0x80 -F 2304 $out`

mapped_pairs=`samtools view -q $qualthresh -f 0x40 -F 2316 $out | eval $mate_filter_cmd | wc -l`
pcr_dupe_pairs=`samtools view -q $qualthresh -u -f 0x40 -F 2316 $out | samtools view -f 0x400 | eval $mate_filter_cmd | wc -l`

mapped_nondupe_pairs=`samtools view -q $qualthresh -f 0x40 -F 3340 $out | eval $mate_filter_cmd | wc -l`
mapped_nondupe_pairs_cis=`samtools view -q $qualthresh -f 0x40 -F 3340 $out | awk '{if (sqrt($9^2) > 0) { print; }}' | eval $mate_filter_cmd | wc -l`
mapped_nondupe_pairs_cis_lt1000=`samtools view -q $qualthresh -f 0x40 -F 3340 $out | awk '{if ((sqrt($9^2) <= 1000) && ($9 != 0)) { print; }}' | eval $mate_filter_cmd | wc -l`
mapped_nondupe_pairs_cis_gt1000=`samtools view -q $qualthresh -f 0x40 -F 3340 $out | awk '{if (sqrt($9^2) > 1000) { print; }}' | eval $mate_filter_cmd | wc -l`
mapped_nondupe_pairs_cis_gt10000=`samtools view -q $qualthresh -f 0x40 -F 3340 $out | awk '{if (sqrt($9^2) > 10000) { print; }}' | eval $mate_filter_cmd | wc -l`
mapped_trans_pairs=`samtools view -q $qualthresh -f 0x40 -F 3340 $out | awk '{if (sqrt($9^2) == 0) { print; }}' | eval $mate_filter_cmd | wc -l`

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

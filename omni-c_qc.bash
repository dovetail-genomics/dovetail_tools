#!/usr/bin/env bash

ref=$1
fq1=$2
fq2=$3
out=$4
rg=$5

bwa mem -5SP -T0 -t4 \
    -R "@RG\tID:$rg\tSM:$rg\tLB:$rg\tPL:ILLUMINA\tPU:none" \
    $ref \
    $fq1 \
    $fq2 \
| samblaster -i stdin -o stdout \
| samtools view -S -h -b \
| samtools sort --threads 4 - > $out ;

samtools index $out ;

preseq lc_extrap -B -P -e 2.1e9 -s 1e8 -seg_len 1000000000 -o $out.preseq $out

ps300m=`cat $out.preseq | grep -P "^300000000.0" | awk '{print $2}'`
 
gt10000=`samtools view -f 64 -F 3328 $out | awk '{if ($9 > 10000) { print; }}' | wc -l`
gt1000=`samtools view -f 64 -F 3328 $out | awk  '{if ($9 > 1000) { print; }}' | wc -l`
gt0=`samtools view -f 64 -F 3328 $out | awk  '{if ($9 > 0) { print; }}' | wc -l`
trans=`samtools view -f 64 -F 3328 $out | awk  '{if ($9 == 0) { print; }}' | wc -l`

echo; echo;
echo "valid pairs          :" `expr $trans + $gt1000`
echo "cis pairs            :" $gt0
echo "cis pairs >1000 bp   :" $gt1000 "(" $(($gt1000 / $gt0)) "% of total cis reads)"
echo "cis pairs >10000 bp  :" $gt10000  "(" $(($gt10000 / $gt0)) "% of total cis reads)"
echo "trans pairs          :" $trans 

echo "Expected unique pairs at 300M sequencing: " $ps300m

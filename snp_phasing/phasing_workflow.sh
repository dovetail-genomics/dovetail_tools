#!/bin/bash

vcf=$1
bam=$2
prefix=$3

#Phase these variants with HAPCUT2

extractHAIRS  --hic 1 --bam ${bam} --VCF ${vcf} --out ${prefix}_hapcut.fragments
HAPCUT2 --hic 1 --fragments ${prefix}_hapcut.fragments --VCF ${vcf} --output ${prefix}_hapcut_output --outvcf 1 

bgzip -c ${prefix}_hapcut_output.phased.VCF  > ${prefix}_hapcut_output.phased.VCF.gz
tabix ${prefix}_hapcut_output.phased.VCF.gz


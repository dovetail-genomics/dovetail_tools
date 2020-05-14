#!/bin/bash

vcf=$1
bam=$2
prefix=$3

#Phase these variants with HAPCUT2

grep -E '^#|0/0|1/1|0/1|1/0|0/2|2/0'  ${vcf} > ${prefix}_filtered_variants.vcf

extractHAIRS  --hic 1 --bam ${bam} --VCF ${prefix}_filtered_variants.vcf --out ${prefix}_hapcut.fragments
HAPCUT2 --hic 1 --fragments ${prefix}_hapcut.fragments --VCF ${prefix}_filtered_variants.vcf --output ${prefix}_hapcut_output --outvcf 1 

bgzip -c ${prefix}_hapcut_output.phased.VCF  > ${prefix}_hapcut_output.phased.VCF.gz
tabix ${prefix}_hapcut_output.phased.VCF.gz


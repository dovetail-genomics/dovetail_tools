#!/bin/bash

if [  $# -le 2 ]
then
        echo "Too few arguments. Please provide all the required arguments."
        echo "Usage: ./phasing_workflow.sh <variants_vcf> <alignment_bam> <output_prefix>"
        exit 1
fi

vcf=$1
bam=$2
prefix=$3

#Phase these variants with HAPCUT2


extractHAIRS  --hic 1 --bam ${bam} --VCF ${vcf} --out ${prefix}_hapcut.fragments
HAPCUT2 --hic 1 --fragments ${prefix}_hapcut.fragments --VCF ${vcf} --output ${prefix}_hapcut_output --outvcf 1 

bgzip -c ${prefix}_hapcut_output.phased.VCF  > ${prefix}_hapcut_output.phased.VCF.gz
tabix ${prefix}_hapcut_output.phased.VCF.gz


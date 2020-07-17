#!/bin/bash


if [  $# -le 5]
then
        echo "Too few arguments. Please provide all the required arguments."
        echo "Usage: ./make_hq_region.sh <reference> <truth_snps.vcf.gz> <called_snps.vcf.gz> <high_confident_regions.bed>\
	       	<output_directory> <output_prefix>"
        exit 1
fi

ref=$1
truth=$2
called=$3
regions=$4
output_dir=$5
output_prefix=$6


gatk Concordance \
        -R ${ref} \
        -isr INTERSECTION \
        -L ${regions} \
        -truth ${truth} \
        -eval ${called} \
        --summary ${output_dir}/${output_prefix}_summary.tsv \
        -tpfp  ${output_dir}/${output_prefix}_tpfp.vcf \
        -tpfn  ${output_dir}/${output_prefix}_tpfn.vcf

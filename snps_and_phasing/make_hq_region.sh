#!/bin/bash

if [  $# -le 4 ]
then
        echo "Too few arguments. Please provide all the required arguments."
        echo "Usage: ./make_hq_region.sh <aligment_bam> <snps_vcf.gz> <GIAB_confident_regions_bed> <output_prefix> <output_directory>"
        exit 1
fi

DOVETOOLS=`dirname $0`

bam=$1
snps=$2
hq_region=$3
output_prefix=$4
output_dir=$5


${DOVETOOLS}/get_HQ_region_bed.py -bam ${bam}  -bedroot ${output_dir}/${output_prefix}
bedtools intersect -a ${hq_region} -b ${output_dir}/${output_prefix}_highconf.bed | grep -v '*' >  ${output_dir}/ConfidentRegionsAND${output_prefix}HQ.bed
bedtools intersect -header -a ${snps} -b ${output_dir}/ConfidentRegionsAND${output_prefix}HQ.bed > ${output_dir}/${output_prefix}_highconf.vcf
bgzip ${output_dir}/${output_prefix}_highconf.vcf
tabix -p vcf ${output_dir}/${output_prefix}_highconf.vcf.gz

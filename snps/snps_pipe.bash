#!/usr/bin/env bash 


if [  $# -le 3 ]
then
        echo "Too few arguments. Please provide all the required arguments."
        echo "Usage: ./snps_pipe.bash <sample_bam> <output_root> <intervals> <reference>"
        exit 1
fi

# Name of sorted indexed bam file.
sample_bam=$1
# Root name/directory for all output files
output_root=$2
# bed file describing the intervals over which snps will be called.
# Typically this would be a confident regions file of some sort, ConfidentRegions.gz
snpcalling_intervals=$3
# e.g. GRCh38.p12.fa
reference = $4

################### Call SNPS #######################

output_vcf=${output_root}.vcf
./hap_caller.bash ${sample_bam} ${output_vcf} ${snpcalling_intervals} ${reference}


################### Filter SNPS #######################

# Create a high confidence SNP file from the intial SNP file
./get_HQ_region_bed.py -bam ${sample_bam} -bedroot ${output_root}
bedtools intersect -a ${output_vcf}.gz -b ${output_root}_highconf.bed > ${output_root}_highconf.vcf 

bgzip ${output_root}_highconf.vcf
tabix -p vcf ${output_root}_highconf.vcf.gz
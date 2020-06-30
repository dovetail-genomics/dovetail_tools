#!/usr/bin/env bash 

if [  $# -le 3 ]
then
        echo "Too few arguments. Please provide all the required arguments."
        echo "Usage: ./hap_caller.bash <sample_bam> <output_vcf> <intervals> <reference>"
        exit 1
fi

# Name of sorted indexed bam file. 
sample_bam=$1
# Name of output vcf file 
output_vcf=$2
# bed file describing the intervals over which snps will be called. 
# Typically this would be a confident regions file of some sort. 
snpcalling_intervals=$3
# e.g. /local/ref/hg38/GRCh38.p12.fa
reference = $4

# omni hc regions Specify the regions of high quality mapping
# observed with this alignment to produce a higher confidence SNP
# output set.   

mkdir -p /local/output/tmp

gatk --java-options -Xmx4096m HaplotypeCaller \
     -I $sample_bam \
     -L $snpcalling_intervals \
     -O ${output_vcf}.g \
     -ERC GVCF \
     -R ${reference} \
     -ip 1000 \
     --native-pair-hmm-threads 1 \
     --max-alternate-alleles 3 \
     -contamination 0 \
     --tmp-dir output/tmp/ \
     --min-base-quality-score 0 

bgzip ${output_vcf}.g
tabix -p vcf ${output_vcf}.g.gz

#     --dbsnp /local/ref/hg38/dbSNP/00-common_all.vcf.gz
gatk GenotypeGVCFs \
     -R ${reference} \
     -V ${output_vcf}.g.gz \
     -L $target_interval
     -O ${output_vcf} -OVI -OVM

bgzip ${output_vcf}
tabix -p vcf ${output_vcf}.gz


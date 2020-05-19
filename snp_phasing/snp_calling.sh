#!/bin/bash


if [  $# -le 2 ]
then
	echo "Too few arguments. Please provide all the required arguments."
	echo "Usage: ./snp_calling.sh <reference_fasta> <alignment_bam> <output_prefix>"
	exit 1
fi


ref=$1
bam=$2
prefix=$3


#Generate contig size file
samtools faidx ${ref}
cut -f 1,2 ${ref}.fai > ${prefix}.contig_size.txt



#Call SNPs with HaplotypeCaller

gatk HaplotypeCaller -A RMSMappingQuality -R ${ref} -I ${bam}  -O ${prefix}_variants.vcf \
	--max-alternate-alleles 3 --contamination-fraction-to-filter 0.0   --minimum-mapping-quality 60 --heterozygosity 0.001 \
	--indel-heterozygosity 1.25E-4 -standard-min-confidence-threshold-for-calling 30.0 --max-genotype-count 1024 --sample-ploidy 2 \
	--num-reference-samples-if-no-call 0 --native-pair-hmm-threads 32


bgzip -c  ${prefix}_variants.vcf > ${prefix}_variants.vcf.gz 
tabix ${prefix}_variants.vcf.gz


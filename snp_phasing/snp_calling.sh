#!/bin/bash

ref=$1
bam=$2
prefix=$3

SRCDIR=`dirname $0`

#Generate contig size file
samtools faidx ${ref}
cut -f 1,2 ${ref}.fai > ${prefix}.contig_size.txt



#Call SNPs with HaplotypeCaller

gatk HaplotypeCaller -A RMSMappingQuality -R ${ref} -I ${bam} -O ${prefix}_raw_variants.vcf -ERC GVCF --native-pair-hmm-threads 1 --max-alternate-alleles 3 --QUIET 

bgzip -c  ${prefix}_variants.vcf > ${prefix}_variants.vcf.gz 
tabix ${prefix}_variants.vcf.gz


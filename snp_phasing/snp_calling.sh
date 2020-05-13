#!/bin/bash

ref=$1
bam=$2
prefix=$3

SRCDIR=`dirname $0`

#Generate contig size file
samtools faidx ${ref}
cut -f 1,2 ${ref}.fai > ${prefix}.contig_size.txt



#Call SNPs with HaplotypeCaller

gatk HaplotypeCaller -R ${ref} -I ${bam} -O ${prefix}_raw_variants.vcf -ERC GVCF --native-pair-hmm-threads 1 --max-alternate-alleles 3 --QUIET 
grep -E '^#|0/0|1/1|0/1|1/0|0/2|2/0'  ${prefix}_raw_variants.vcf > ${prefix}_variants.vcf

#First call SNPs with UnifiedGenotyper

#java -jar ${SRCDIR}/GenomeAnalysisTK.jar  -T UnifiedGenotyper -R ${ref} -drf BadMate \
#       -I ${bam} -o ${prefix}_variants.vcf -L ${confident_regions}


bgzip -c  ${prefix}_variants.vcf > ${prefix}_variants.vcf.gz 
tabix ${prefix}_variants.vcf.gz

python ${SRCDIR}/get_HQ_region.py -bam ${bam} -qthresh 0   -bedroot ${prefix}  -contigsizes ${prefix}.contig_size.txt

bedtools merge -i ${prefix}.bed > ${prefix}_merged.bed
bedtools complement -i ${prefix}_merged.bed -g ${prefix}.contig_size.txt  > ${prefix}_highconf.bed


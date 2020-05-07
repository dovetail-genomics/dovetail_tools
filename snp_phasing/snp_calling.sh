#!/bin/bash

ref=$1
bam=$2
confident_regions=$3
prefix=$4

SRCDIR=`dirname $0`

#Generate contig size file
samtools faidx ${ref}
cut -f 1,2 ${ref}.fai > ${prefix}.contig_size.txt


#First call SNPs with UnifiedGenotyper

java -jar ${SRCDIR}/GenomeAnalysisTK.jar  -T UnifiedGenotyper -R ${ref} -drf BadMate \
       -I ${bam} -o ${prefix}_variants.vcf -L ${confident_regions}

bgzip -c ${prefix}_variants.vcf > ${prefix}_variants.vcf.gz 
tabix ${prefix}_variants.vcf.gz
#Now filter these SNPs using custon scripts

bgzip -c  ${prefix}_variants.vcf > ${prefix}_variants.vcf.gz 
tabix ${prefix}_variants.vcf.gz

python ${SRCDIR}/qd_fs_filter.py -vcfin ${prefix}_variants.vcf.gz -qdthresh 2 -fsthresh 60 -vcfout ${prefix}_variants_qd2_fs60.vcf


python ${SRCDIR}/get_HQ_region.py -bam ${bam} -qthresh 60   -bedroot ${prefix}  -contigsizes ${prefix}.contig_size.txt

bedtools merge -i ${prefix}.bed > ${prefix}_merged.bed
bedtools complement -i ${prefix}_merged.bed -g ${prefix}.contig_size.txt  > ${prefix}_highconf.bed

bedtools intersect -a ${prefix}_variants_qd2_fs60.vcf.gz  -b ${prefix}_highconf.bed -header > ${prefix}_variants_qd2_fs60_highconf.vcf

bgzip -c ${prefix}_variants_qd2_fs60_highconf.vcf > ${prefix}_variants_qd2_fs60_highconf.vcf.gz 
tabix ${prefix}_variants_qd2_fs60_highconf.vcf.gz

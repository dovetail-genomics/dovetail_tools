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

java -jar GenomeAnalysisTK.jar  -T UnifiedGenotyper -R ${ref} -drf BadMate \
       -I ${bam} -o ${prefix}_variants.vcf -L ${confident_regions}


#Now filter these SNPs using custon scripts

python ${SRCDIR}/qd_fs_filter.py -vcfin ${prefix}_variants.vcf -qdthresh 2 -fsthresh 60 -vcfout ${prefix}_variants_qd2_fs60.vcf

python ${SRCDIR}/get_LQ_region_bed.py -bam ${bam} -qthresh 60   -bedroot ${prefix}  -contigsizes ${prefix}.contig_size.txt

bedtools intersect -a ${prefix}_variants_qd2_fs60.vcf -b ${prefix}_highconf.bed -header > ${prefix}_variants_qd2_fs60_highconf.vcf

bgzip -c ${prefix}_variants_qd2_fs60_highconf.vcf > ${prefix}_variants_qd2_fs60_highconf.vcf.gz
tabix ${prefix}_variants_qd2_fs60_highconf.vcf.gz


#Phase these variants with HAPCUT2

extractHAIRS  --hic 1 --bam ${bam} --VCF ${prefix}_variants_qd2_fs60_highconf.vcf --out ${prefix}_hapcut.fragments
HAPCUT2 --hic 1 --fragments ${prefix}_hapcut.fragments --VCF ${prefix}_variants_qd2_fs60_highconf.vcf --output hapcut_output --outvcf 1 

bgzip -c hapcut_output.phased.VCF  > hapcut_output.phased.VCF.gz
tabix hapcut_output.phased.VCF.gz


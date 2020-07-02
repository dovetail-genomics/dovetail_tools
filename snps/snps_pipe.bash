#!/usr/bin/env bash 

# NOTE:  Assumes you have put dovetail_tools/snps/ in your path. 
# e.g. 
# export PATH=$PATH:/local/ubuntu/src/dovetail_tools/snps/
#

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
reference=$4

################### Call SNPS #######################

output_vcf=${output_root}.vcf
call_snps.bash ${sample_bam} ${output_vcf} ${snpcalling_intervals} ${reference}


################### Filter SNPS #######################

# Create a high confidence SNP file from the intial SNP file
get_HQ_region_bed.py -bam ${sample_bam} -bedroot ${output_root} 
bedtools intersect -header -a ${output_vcf}.gz -b ${output_root}_highconf.bed > ${output_root}_highconf.vcf 
bgzip ${output_root}_highconf.vcf
tabix -p vcf ${output_root}_highconf.vcf.gz



# Note:  If you want to compute SNP concordance with gatk Concordance you will
# need to create a high confident bed file that is the intersection of the ConfidentRegions.bed
# that ships with, for example, Illumina Platinum truthset AND the OmniC confident regions
# bed file created above, ${output_root}_highconf.bed and pass that into gatkConcordance as 
# the interval like: 

# bedtools intersect -a ${snpcalling_intervals} -b ${output_root}_highconf.bed > ConfidentANDOmniConfident.bed
#gatk Concordance \
#     -R /local/ref/hg38/GRCh38.p12.fa \
#     -isr INTERSECTION \
#     -L  ConfidentANDOmniConfident.bed  \
#     -truth $truthset \
#     -eval $eval_vcf \
#     --summary ${output_root}_summary.tsv \
#     -tpfp ${output_root}_tpfp.vcf \
#     -tpfn ${output_root}_tpfn.vcf

# This will restrict the evaltuion to the confident regions and give you a clearer picture of how good
# the confident SNP calls probably are. 

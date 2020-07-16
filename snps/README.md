
# Calling SNPs with OmniC data

SNPs can be called by editing the makefile to point to a base directory that contains the snp resources and a directory that contains dovetail_tools

'''
ROOT=/home/ubuntu/prj/snpexample/
DOVETOOLS=/home/ubuntu/src/dovetail_tools/
LIB_NAME=LM979
'''

A .json file is expected in:

'''
${ROOT}/${LIB_NAME}.json
'''

which points to the bam file, like so:

'''
{
    "germline":[
        {
            "bams":[
             "/home/ubuntu/prj/snpexample/bam/LM979.bam"
            ]
        }
    ]
}
'''

The directory, here snpexample, will need to contain the following genome reference files (or symlinks to them):

'''
snpexample/ref/hg38/GRCh38.p12.fa
snpexample/ref/hg38/GRCh38.p12.fa.fai
snpexample/ref/hg38/GRCh38.p12.fa.sizes
snpexample/ref/hg38/dbSNP/00-common_all.vcf.gz
snpexample/ref/hg38/dbSNP/00-common_all.vcf.gz.tbi
'''

Where the fasta file is the reference fasta file and the vcfs under dbSNP are the commonly occurring SNPs from dbSNP.  


This snp pipeline  produces a raw set of snps and a separate high confidence set of snps.  High confident snps are snps that fall in regions where there are no MQ=0 aligned reads.  Some regions that are confident regions for shotgun are not confident regions for OmniC because the pair size distribution of *-C data is very broad and some mates whose placements could be disambiguated by bwa with shotgun data can not be disambiguated with *-C data.   A modified version of bwa, or some post-processing, might be able to constrain some of these mates based on the broader Omni-C pair size distribution but for now we just identify these difficult regions and exclude them from the results. 

  







# Calling SNPs with OmniC data

Before you run SNP analysis, make sure you have conda environment set up from the yaml file in the `conda` folder. To call SNPs, you can run `snp_pipeline.py` script as follows:

```
usage: snps_pipeline.py [-h] -B BAM -O OUT_DIR -R REF -L INTERVALS -dbsnp
                        DBSNP [-clean] [-t MAX_CORES] [-drm DRM]

optional arguments:
  -h, --help            show this help message and exit
  -B BAM, --bam BAM
  -O OUT_DIR, --out_dir OUT_DIR
  -R REF, --ref REF
  -L INTERVALS, --intervals INTERVALS
                        Regions to call intervals on in bed format. Helpful for
                        parallelization as Each interval is processed
                        independently
  -dbsnp DBSNP, --dbsnp DBSNP
  -clean, --clean       Removes split vcf files
  -t MAX_CORES, --max_cores MAX_CORES
  -drm DRM
```

`-L` option is used for multiprocessing by calling SNPs on each interval independently. You can use `hg38_intervals.txt` file provided in the repository to use with this option. `-dbsnp` requires VCF file for dbSNP database SNPs and it's index. This is required for Base Quality Recalibrator and Genotyping SNPs. You can obtain there here: https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/GATK/. `-t` option can be used to specify number of cores. With 8 cores, the snp pipeline takes about 24 hours to generate final VCF file. 

This snp pipeline  produces a raw set of snps and a separate high confidence set of snps.  High confident snps are snps that fall in regions where there are no MQ=0 aligned reads.  Some regions that are confident regions for shotgun are not confident regions for OmniC because the pair size distribution of *-C data is very broad and some mates whose placements could be disambiguated by bwa with shotgun data can not be disambiguated with *-C data.   A modified version of bwa, or some post-processing, might be able to constrain some of these mates based on the broader Omni-C pair size distribution but for now we just identify these difficult regions and exclude them from the results. We can generate these regions from `make_hq_regions.sh` script as follows:

```
Usage: ./make_hq_region.sh <aligment_bam> <snps_vcf.gz> <GIAB_confident_regions_bed> <output_prefix> <output_directory>
```
The GIAB confident region bed file for shotgun data can be obtained from here: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/. The output will be a bed file describing the confident regions along with the filtered VCF file based on these regions. 

To evaluate SNPs against the truthset, you can use `snp_concordance.sh` script as follows:
```
./make_hq_region.sh <reference> <truth_snps.vcf.gz> <called_snps.vcf.gz> <high_confident_regions.bed> <output_directory> <output_prefix>
```

Here, for the `<high_confident_regions.bed>` option, use the bed file obtained from the `make-hq_region.sh` script and use the filtered VCF file for `<called_snps.vcf.gz>` option. 


<!---
SNPs can be called by editing the makefile to point to a base directory that contains the snp resources and a directory that contains dovetail_tools

```
ROOT=/home/ubuntu/prj/snpexample/
DOVETOOLS=/home/ubuntu/src/dovetail_tools/
LIB_NAME=LM979
```

The directory, here snpexample, will need to contain the following genome reference files (or symlinks to them):

```
snpexample/ref/hg38/GRCh38.p12.fa
snpexample/ref/hg38/GRCh38.p12.fa.fai
snpexample/ref/hg38/GRCh38.p12.fa.sizes
snpexample/ref/hg38/dbSNP/00-common_all.vcf.gz
snpexample/ref/hg38/dbSNP/00-common_all.vcf.gz.tbi
```

Where the fasta file is the reference fasta file and the vcfs under dbSNP are the commonly occurring SNPs from dbSNP.  

To filter low quality regions a ConfidentRegions.bed.gz is expected in:

```
snpexample/truthset/ConfidentRegions.bed.gz
snpexample/truthset/NA12878.vcf.gz
```

The locations and names of these files can be altered in the included makefile. The max cpus allowed can also be edited in the makefile. To call snps, set up these files and then do:

```
make all
```

Or to run each step separately:

```
make snp_pipe
make hqregion
make concordance
```







This snp pipeline  produces a raw set of snps and a separate high confidence set of snps.  High confident snps are snps that fall in regions where there are no MQ=0 aligned reads.  Some regions that are confident regions for shotgun are not confident regions for OmniC because the pair size distribution of *-C data is very broad and some mates whose placements could be disambiguated by bwa with shotgun data can not be disambiguated with *-C data.   A modified version of bwa, or some post-processing, might be able to constrain some of these mates based on the broader Omni-C pair size distribution but for now we just identify these difficult regions and exclude them from the results. 

 -->
# Phasing SNPs with Omni-C data 
## Environment

In order to run these scripts, you would require following dependencies:

- [tabix](https://anaconda.org/bioconda/tabix)
- [WhatsHap](https://whatshap.readthedocs.io/en/latest/)
- [HapCUT2](https://github.com/vibansal/HapCUT2)

All the dependencies can be installed with Conda. Also, we assume you have all the depencies required for alignment/qc pipeline installed at this time. 

## Phasing variants with Omni-C data

To phase variants, we use HapCUT2 phasing program. The wrapper script to run for phasing is `phasing_workflow.sh`. It can be run as

```
./phasing_workflow.sh <variants_vcf> <alignment_bam> <output_prefix>
```


Note that `<vcf>` file here should NOT be gzipped as HapCUT2 requires unzipped VCF file. This script will first filter input VCF to contain just signel or biallelic variants as HapCUT2 operataes only on these types of variants. The output of the script is `<output_prefix>_hapcut_output.phased.VCF.gz` file that contains phased variants. 

The phasing can be evaluated if the ground truth is known. This can be done with `phasing_eval.sh` script. You can do this as follow: 

```
./phasing_eval.sh <output_vcf.gz> <truth_vcf.gz> <output_prefix>
```


# Calling SNPs with OmniC data

Snps can be called with the included scripts like: 

```
snps_pipe.bash sample.bam outdir/rootname ConfidentRegions.bed.gz GRCh38.p12.fa 
```

Where 'rootname' is the root name for all the output files and ConfidentRegions.bed is a file describing the intervals over which SNPs should be called.  

This produces a raw set of snps and a separate high confidence set of snps.  High confident snps are snps that fall in regions where there are no MQ=0 aligned reads.  Some regions that are confident regions for shotgun are not confident regions for OmniC because the pair size distribution of *-C data is very broad and some mates whose placements could be disambiguated by bwa with shotgun data can not be disambiguated with *-C data.   A modified version of bwa, or some post-processing, might be able to constrain some of these mates based on the broader Omni-C pair size distribution but for now we just identify these difficult regions and exclude them from the results. 

 
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

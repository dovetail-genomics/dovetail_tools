# Omni-C SNV calling and phasing pipeline

In order to run this pipeline, you need to align and QC data using `omni-c_qc.bash` script from this repository. Once you have the BAM file, you can proceed with the SNP calling. 

## Environment

In order to run these scripts, you would require following dependencies:

- [tabix](https://anaconda.org/bioconda/tabix)
- [WhatsHap](https://whatshap.readthedocs.io/en/latest/)
- [HapCUT2](https://github.com/vibansal/HapCUT2)
- [java](https://anaconda.org/bioconda/java-jdk)
- [pysam](https://pysam.readthedocs.io/en/latest/index.html)

All the dependencies can be installed with Conda. Also, we assume you have all the depencies required for alignment/qc pipeline installed at this time. 

## SNV calling with Omni-C data

We use `HaplotypeCaller` program from GATK to call SNVs.

```
./snp_calling.sh <reference_fasta> <alignment_bam> <output_prefix>
```
This run will produce `vcf`, `vcf.gz`, and `vcf.gz.tbi` file in the end. 

## Phasing variants with Omni-C data

To phase variants, we use HapCUT2 phasing program. The wrapper script to run for phasing is `phasing_workflow.sh`. It can be run as

```
./phasing_workflow.sh <variants_vcf> <alignment_bam> <output_prefix>
```


Note that `<vcf>` file here should NOT be gzipped as HapCUT2 requires unzipped VCF file. This script will first filter input VCF to contain just signel or biallelic variants as HapCUT2 operataes only on these types of variants. The output of the script is `<output_prefix>_hapcut_output.phased.VCF.gz` file that contains phased variants. 

## Evaluating phasing

If you want to evaluate phasing against ground truth/trio phased variant calls, you can use `phasing_eval.sh` script. It uses WhatsHap program to generate summary statistics. You can run it as

```
./phasing_eval <hapcut_phasing_vcf> <truth_phasing_vcf> <output_prefix>
```

Note that vcf files here should be gzipped and indexed. 

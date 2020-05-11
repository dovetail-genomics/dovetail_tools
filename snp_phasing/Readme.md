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

We use `UnifiedGenotyper` program from GATK to call SNVs. Since it is not available in the latest build of GATK, we provide the exact GATK jar file in this folder that contins `UnifiedGenotyper`. The wrapper script that you need to run to call SNVs is `snp_calling.sh` and you need to run it as follows:

```
./snp_calling.sh <reference_fasta_file> <bam_file> <confident_region_bed_file> <output_prefix>
```

Here, the `confident_region_bed_file` is the bed file of high confidence regions released by GIAB consortium. We also do additional filtering of variants based on strand bias. Furthermore, we also ignore variants in the regions where reads map with mapping quality 0. `snp_calling.sh` takes care of all the filtering and it would produce `vcf`, `vcf.gz`, and `vcf.gz.tbi` file in the end. 

If you want to evaluate variant calling accuracy with some groundtruth of known variants, you can use `snp_eval.sh` script. It runs `GenotypeConcordance` from GATK. It can be run as:

```
./snp_eval.sh <reference_fasta_file> <truthset_variants> <called_variants> <output_prefix>
```

Note that `<truthset_variants>` and `<called_variants>` should be gzipped and indexed. 

## Phasing variants with Omni-C data

To phase variants, we use HapCUT2 phasing program. The wrapper script to run for phasing is `phasing_workflow.sh`. It can be run as

```
./phasing_workflow.sh <vcf> <bam> <output_prefix>
```

Note that `<vcf>` file here should NOT be gzipped as HapCUT2 requires unzipped VCF file. The output of the script is `<output_prefix>_hapcut_output.phased.VCF.gz` file that contains phased variants. 

## Evaluating phasing

If you want to evaluate phasing against ground truth/trio phased variant calls, you can use `phasing_eval.sh` script. It uses WhatsHap program to generate summary statistics. You can run it as

```
./phasing_eval <hapcut_phasing_vcf> <truth_phasing_vcf> <output_prefix>
```

Note that vcf files here should be gzipped and indexed. 
